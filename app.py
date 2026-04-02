"""
TB Dashboard - A Dash application for visualizing TB genomic data and drug resistance.
"""
import dash
from dash import html, dcc, Input, Output, State, callback_context, dash_table
import dash_bootstrap_components as dbc
import dash_jbrowse
import pandas as pd
from typing import Optional
import os

from data_utils import DataLoader
from coordinate_calculator import CoordinateCalculator

# Initialize data loader
data_loader = DataLoader(data_dir="data")
data_loader.load_gff3()
data_loader.load_catalogue()
data_loader.load_genomic_coordinates()

coord_calculator = CoordinateCalculator(data_loader)

from flask import send_from_directory

# Initialize Dash app
app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.BOOTSTRAP, dbc.icons.BOOTSTRAP],
    suppress_callback_exceptions=True
)
app.title = "TB Dashboard - Genomic Resistance Explorer"

# Flask route to serve data files for JBrowse
@app.server.route('/data/<path:path>')
def serve_data(path):
    return send_from_directory('data', path)

# App layout
app.layout = dbc.Container([
    dcc.Store(id="current-gene-store"),
    # Header
    dbc.Row([
        dbc.Col([
            html.H1(
                [dbc.Badge("TB", color="danger", className="me-2"), 
                 "Dashboard - Genomic Resistance Explorer"],
                className="my-4"
            ),
            html.P(
                "Search for genes or mutations to visualize genomic regions and explore drug resistance profiles",
                className="lead text-muted"
            )
        ])
    ], className="mb-4"),
    
    # Search Section
    dbc.Card([
        dbc.CardBody([
            dbc.Row([
                dbc.Col([
                    dbc.Label("Search Gene or Mutation", className="fw-bold"),
                    dbc.InputGroup([
                        dbc.Input(
                            id="search-input",
                            type="text",
                            placeholder="Enter gene name (e.g., dnaA, gyrA, rpoB) or mutation (e.g., dnaA_c.102G>A)",
                            debounce=True
                        ),
                        dbc.Button(
                            html.I(className="bi bi-search"),
                            id="search-button",
                            color="primary",
                            className="ms-2"
                        )
                    ])
                ], width=8),
                dbc.Col([
                    dbc.Label("Quick Searches", className="fw-bold"),
                    dbc.ButtonGroup([
                        dbc.Button("rpoB", id="quick-rpoB", size="sm", className="me-1"),
                        dbc.Button("gyrA", id="quick-gyrA", size="sm", className="me-1"),
                        dbc.Button("katG", id="quick-katG", size="sm"),
                    ])
                ], width=4)
            ])
        ])
    ], className="mb-4"),
    
    # Results Section
    html.Div(id="search-results-container", children=[]),
    
    # Footer
    dbc.Row([
        dbc.Col([
            html.Hr(),
            html.P([
                "TB Dashboard v1.0 | Data sources: H37Rv reference genome, WHO catalogue of mutations"
            ], className="text-muted text-center small")
        ])
    ])
    
], fluid=True, className="py-4")


# Callbacks
@app.callback(
    [Output("search-results-container", "children"),
     Output("current-gene-store", "data")],
    [Input("search-button", "n_clicks"),
     Input("quick-rpoB", "n_clicks"),
     Input("quick-gyrA", "n_clicks"),
     Input("quick-katG", "n_clicks")],
    State("search-input", "value")
)
def search_gene(
    search_clicks: Optional[int],
    rpoB_clicks: Optional[int],
    gyrA_clicks: Optional[int],
    katG_clicks: Optional[int],
    search_value: Optional[str]
):
    """Handle search input and display results."""
    ctx = callback_context
    if not ctx.triggered:
        return [], None
    
    # Determine which button was clicked
    triggered_id = ctx.triggered[0]["prop_id"].split(".")[0]
    
    if triggered_id == "quick-rpoB":
        search_value = "rpoB"
    elif triggered_id == "quick-gyrA":
        search_value = "gyrA"
    elif triggered_id == "quick-katG":
        search_value = "katG"
    
    if not search_value:
        return [], None
    
    # Try to find gene
    gene_info = data_loader.get_gene_info(search_value)
    
    if not gene_info:
        # Try extracting gene from mutation format (e.g., dnaA_c.102G>A)
        if "_" in search_value:
            gene_part = search_value.split("_")[0]
            gene_info = data_loader.get_gene_info(gene_part)
            
    if not gene_info:
        # Try searching for partial matches
        genes = data_loader.search_genes(search_value)
        if len(genes) > 0:
            gene_info = genes[0]
        else:
            return dbc.Alert(f"Gene '{search_value}' not found. Please try another search.", color="warning"), None
    
    # Get gene mutations
    mutations_df = data_loader.search_mutations_by_gene(gene_info.gene_name)
    
    # Get drug resistance info
    drug_resistance = data_loader.get_drug_resistance_info(gene_info.gene_name)
    
    # Store gene_info as dict
    gene_data = {
        "gene_name": gene_info.gene_name,
        "locus_tag": gene_info.locus_tag,
        "start": gene_info.start,
        "end": gene_info.end,
        "strand": gene_info.strand,
        "chromosome": gene_info.chromosome,
        "product": gene_info.product
    }
    
    # Build results
    results = []
    
    # Gene Info Card
    gene_length = gene_info.end - gene_info.start + 1
    results.append(
        dbc.Card([
            dbc.CardHeader([
                dbc.Row([
                    dbc.Col([
                        html.H4(f"{gene_info.gene_name}", className="mb-0"),
                        html.Small(f"Locus: {gene_info.locus_tag}", className="text-muted")
                    ]),
                    dbc.Col([
                        html.Div(f"Product: {gene_info.product or 'Unknown'}", className="text-muted")
                    ])
                ])
            ], className="bg-light"),
            dbc.CardBody([
                dbc.Row([
                    dbc.Col([
                        html.Strong("Genomic Coordinates: "),
                        html.Span(f"{gene_info.chromosome}:{gene_info.start:,}-{gene_info.end:,}")
                    ]),
                    dbc.Col([
                        html.Strong("Strand: "),
                        html.Span("+" if gene_info.strand == "+" else "-")
                    ]),
                    dbc.Col([
                        html.Strong("Gene Length: "),
                        html.Span(f"{gene_length:,} bp")
                    ])
                ])
            ])
        ], className="mb-4")
    )
    
    # JBrowse Card
    jbrowse_start = max(1, gene_info.start - 500)
    jbrowse_end = gene_info.end + 500
    jb_config = data_loader.get_jbrowse_config(gene_info.chromosome, jbrowse_start, jbrowse_end)
    
    results.append(
        dbc.Card([
            dbc.CardHeader([
                html.H5([html.I(className="bi bi-eye"), " Genomic Visualization"], className="mb-0")
            ]),
            dbc.CardBody([
                dash_jbrowse.LinearGenomeView(
                    id="jbrowse-linear-view",
                    assembly=jb_config["assembly"],
                    tracks=jb_config["tracks"],
                    defaultSession=jb_config["defaultSession"],
                    location=f"{gene_info.chromosome}:{jbrowse_start}-{jbrowse_end}"
                ),
                html.Div([
                    dbc.Button(
                        "Open in Full JBrowse",
                        id="open-jbrowse",
                        color="link",
                        className="mt-2",
                        href=f"#jbrowse?region={gene_info.chromosome}:{jbrowse_start}-{jbrowse_end}"
                    )
                ])
            ])
        ], className="mb-4")
    )
    
    # Drug Resistance Card
    if len(drug_resistance) > 0:
        drugs = drug_resistance['drug'].unique()
        drug_cards = []
        for drug in drugs:
            drug_data = drug_resistance[drug_resistance['drug'] == drug]
            tiers = drug_data['tier'].unique()
            
            drug_cards.append(
                dbc.Card([
                    dbc.CardHeader([
                        html.H6(drug, className="mb-0"),
                        dbc.Badge(f"{len(drug_data)} mutations", color="secondary", className="ms-2")
                    ]),
                    dbc.CardBody([
                        html.P([
                            html.Strong("Tiers: "),
                            html.Span(", ".join([str(t) for t in tiers if pd.notna(t)]))
                        ], className="mb-2"),
                        html.Div(
                            id=f"drug-detail-{drug}",
                            children=create_drug_detail_table(drug_data, gene_info.gene_name)
                        )
                    ])
                ], className="mb-2")
            )
        
        results.append(
            dbc.Card([
                dbc.CardHeader([
                    html.H5([html.I(className="bi bi-capsule"), " Drug Resistance Profile"], className="mb-0")
                ], className="bg-light"),
                dbc.CardBody(drug_cards)
            ], className="mb-4")
        )
    
    # Genomic Coordinates Card
    if len(mutations_df) > 0:
        results.append(
            dbc.Card([
                dbc.CardHeader([
                    html.H5([html.I(className="bi bi-map"), " Genomic Coordinates Details"], className="mb-0"),
                    html.Small(f"{len(mutations_df)} unique mutations found", className="text-muted ms-2")
                ], className="bg-light"),
                dbc.CardBody([
                    create_genomic_coords_table(mutations_df, gene_info)
                ])
            ], className="mb-4")
        )
    
    # Calculator Card (placeholder for coordinate calculations)
    results.append(
        dbc.Card([
            dbc.CardHeader([
                html.H5([html.I(className="bi bi-calculator"), " Coordinate Calculator"], className="mb-0")
            ], className="bg-info bg-opacity-10"),
            dbc.CardBody([
                html.P("Select a mutation from the table above to see detailed coordinate calculations."),
                html.Div(id="coordinate-calculator-display")
            ])
        ], className="mb-4")
    )
    
    return results, gene_data


@app.callback(
    Output("coordinate-calculator-display", "children"),
    [Input("genomic-coords-table", "active_cell")],
    [State("genomic-coords-table", "derived_viewport_data"),
     State("current-gene-store", "data")],
    prevent_initial_call=True
)
def display_coordinate_calculation(active_cell, viewport_data, gene_data):
    """Update coordinate calculator display based on selected mutation."""
    if not active_cell or not viewport_data or not gene_data:
        return html.P("Select a mutation from the table above to see detailed coordinate calculations.")
    
    # Get the row data from the current viewport (handles pagination/filtering)
    row_idx = active_cell['row']
    
    if row_idx >= len(viewport_data):
        return html.P("Error: Selected row not found in viewport.")
        
    row = viewport_data[row_idx]
    variant = row['variant']
    gene_name = gene_data['gene_name']
    
    # Use the SPECIFIC genomic position from the clicked row
    actual_genomic_pos = int(row['position']) if str(row['position']).isdigit() else None
    ref = row.get('reference_nucleotide', 'N/A')
    alt = row.get('alternative_nucleotide', 'N/A')
    
    # Calculate relative details for display context
    # Note: coord_calculator.calculate_relative_position is now in CoordinateCalculator
    from data_utils import GeneInfo
    rel_pos = coord_calculator.calculate_relative_position(actual_genomic_pos, 
                                                         GeneInfo(
                                                             gene_id="", gene_name=gene_name, locus_tag=gene_data['locus_tag'],
                                                             start=gene_data['start'], end=gene_data['end'], strand=gene_data['strand'],
                                                             product=gene_data['product'], chromosome=gene_data['chromosome']
                                                         )) if actual_genomic_pos else None
    
    # Get ALL drug resistance entries for this variant
    resistance_df = data_loader.get_drug_resistance_info(gene_name, variant)
    
    # If no results found for specific fs/ins/del, try LoF if it's a loss-of-function type
    if len(resistance_df) == 0:
        if any(term in variant.lower() for term in ["fs", "stop", "ins", "del"]):
            resistance_df = data_loader.get_drug_resistance_info(gene_name, "LoF")
    
    gen_pos_display = f"{actual_genomic_pos:,}" if actual_genomic_pos else "N/A"
    
    # Create Drug Resistance Summary
    if len(resistance_df) > 0:
        res_items = []
        for _, res_row in resistance_df.iterrows():
            res_items.append(html.Div([
                dbc.Row([
                    dbc.Col(html.Span([
                        html.Strong(res_row['drug']),
                        html.Br(),
                        html.Small(f"({res_row['gene']})", className="text-muted")
                    ]), width=4),
                    dbc.Col(dbc.Badge(res_row.get('FINAL CONFIDENCE GRADING', 'N/A'), 
                                     color="danger" if "Assoc w R" in str(res_row.get('FINAL CONFIDENCE GRADING', '')) else "secondary"), 
                           width=5),
                    dbc.Col(html.Small(f"T{res_row.get('tier', 'N/A')}"), width=3),
                ], className="py-1 border-bottom align-items-center")
            ]))
        res_display = html.Div(res_items, className="mt-2")
    else:
        res_display = html.P("No drug resistance data found in catalogue for this variant.", className="text-muted small")

    # Create display
    return html.Div([
        dbc.Row([
            dbc.Col([
                html.H5([html.I(className="bi bi-info-circle me-2"), variant], className="text-primary"),
                html.Hr()
            ])
        ]),
        dbc.Row([
            dbc.Col([
                html.Strong("Genomic Position: "), html.Span(gen_pos_display),
                html.Br(),
                html.Strong("Nucleotide Change: "), html.Span(f"{ref} > {alt}"),
                html.Br(),
                html.Strong("Gene Relative (c.): "), html.Span(f"c.{rel_pos}" if rel_pos else "N/A")
            ], width=6),
            dbc.Col([
                html.Strong("Gene: "), html.Span(f"{gene_name} ({gene_data['strand']})"),
                html.Br(),
                html.Strong("Coordinates: "), html.Small(f"{gene_data['start']:,}-{gene_data['end']:,}")
            ], width=6)
        ], className="mb-3"),
        
        dbc.Row([
            dbc.Col([
                html.H6("Drug Resistance Profile", className="border-bottom pb-1"),
                res_display
            ])
        ], className="mb-3"),
        
        html.Div([
            dbc.Accordion([
                dbc.AccordionItem([
                    html.Code(
                        f"""
Strand: {gene_data['strand']}
Gene Start: {gene_data['start']:,}
Gene End: {gene_data['end']:,}

Formula ({gene_data['strand']} strand):
c.pos = genomic - start + 1 (for +)
c.pos = end - genomic + 1 (for -)

Calculation:
{f"{actual_genomic_pos:,} - {gene_data['start']:,} + 1" if gene_data['strand'] == '+' else f"{gene_data['end']:,} - {actual_genomic_pos:,} + 1"} = {rel_pos}
                        """,
                        style={"whiteSpace": "pre-wrap", "display": "block", "padding": "10px", "backgroundColor": "#f8f9fa"}
                    )
                ], title="Coordinate Calculation Details")
            ], start_collapsed=True)
        ])
    ])


def create_drug_detail_table(drug_data: pd.DataFrame, gene_name: str) -> html.Div:
    """Create detailed drug resistance table."""
    # Show top mutations by tier
    tier_order = {'1': 1, '2': 2, '3': 3, '4': 4, '5': 5}
    drug_data_sorted = drug_data.copy()
    drug_data_sorted['tier_sort'] = drug_data_sorted['tier'].map(tier_order).fillna(99)
    drug_data_sorted = drug_data_sorted.sort_values('tier_sort')
    
    # Columns we want to display
    cols = [
        {'name': 'Mutation', 'id': 'mutation'},
        {'name': 'Tier', 'id': 'tier'},
        {'name': 'Effect', 'id': 'effect'},
        {'name': 'Confidence', 'id': 'FINAL CONFIDENCE GRADING'}
    ]
    
    table = dash_table.DataTable(
        data=drug_data_sorted.to_dict('records'),
        columns=cols,
        filter_action="native",
        sort_action="native",
        page_size=10,
        style_table={'overflowX': 'auto'},
        style_cell={
            'textAlign': 'left',
            'minWidth': '100px', 'maxWidth': '300px',
            'overflow': 'hidden',
            'textOverflow': 'ellipsis',
            'fontFamily': 'inherit'
        },
        style_header={
            'backgroundColor': '#f8f9fa',
            'fontWeight': 'bold'
        }
    )
    
    return html.Div([
        table,
        html.Small(f"Showing {len(drug_data)} mutations", className="text-muted")
    ])


def create_genomic_coords_table(mutations_df: pd.DataFrame, gene_info) -> html.Div:
    """Create genomic coordinates drill-down table."""
    display_df = mutations_df.copy()
    
    # Add calculated columns
    display_df['Gene Relative Pos'] = display_df.apply(
        lambda row: coord_calculator.calculate_relative_position(
            int(row['position']), gene_info
        ) if row.get('position') and str(row['position']).isdigit() else None,
        axis=1
    )
    
    # Format gene relative
    display_df['Gene Relative'] = display_df['Gene Relative Pos'].apply(
        lambda x: f"c.{int(x)}" if pd.notna(x) else 'N/A'
    )
    
    # Format position
    display_df['Position Display'] = display_df['position'].apply(
        lambda x: f"{int(x):,}" if str(x).isdigit() else 'N/A'
    )
    
    cols = [
        {'name': 'Variant', 'id': 'variant'},
        {'name': 'Genomic Position', 'id': 'Position Display'},
        {'name': 'Ref', 'id': 'reference_nucleotide'},
        {'name': 'Alt', 'id': 'alternative_nucleotide'},
        {'name': 'Gene Relative', 'id': 'Gene Relative'}
    ]
    
    table = dash_table.DataTable(
        id="genomic-coords-table",
        data=display_df.to_dict('records'),
        columns=cols,
        # Allow row selection by clicking cells
        cell_selectable=True,
        # Synchronize with callback
        filter_action="native",
        sort_action="native",
        page_size=20,
        style_table={'overflowX': 'auto'},
        style_cell={
            'textAlign': 'left',
            'minWidth': '80px', 'maxWidth': '250px',
            'overflow': 'hidden',
            'textOverflow': 'ellipsis',
            'fontFamily': 'inherit'
        },
        style_header={
            'backgroundColor': '#f8f9fa',
            'fontWeight': 'bold'
        },
        # Tooltip for long nucleotide sequences
        tooltip_data=[
            {
                column: {'value': str(value), 'type': 'markdown'}
                for column, value in row.items()
            } for row in display_df.to_dict('records')
        ],
        tooltip_duration=None
    )
    
    return html.Div([
        table,
        html.Small(f"Total entries: {len(mutations_df)}", className="text-muted")
    ])


if __name__ == "__main__":
    print("Starting TB Dashboard...")
    print("Access the application at: http://localhost:8050")
    app.run(debug=True, host="0.0.0.0", port=8050)
