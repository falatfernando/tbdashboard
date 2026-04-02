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

# Get absolute path for data directory
DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

# Initialize data loader
data_loader = DataLoader(data_dir=DATA_DIR)
data_loader.load_gff3()
data_loader.load_catalogue()
data_loader.load_genomic_coordinates()

coord_calculator = CoordinateCalculator(data_loader)

from flask import send_from_directory

# Initialize Dash app
app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.BOOTSTRAP, dbc.icons.BOOTSTRAP],
    suppress_callback_exceptions=True,
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}]
)
app.title = "TB Dashboard | LaPAM"

# Expose the Flask server for gunicorn
server = app.server

# Flask route to serve data files for JBrowse
@server.route('/data/<path:path>')
def serve_data(path):
    return send_from_directory(DATA_DIR, path)

# Navigation Bar
navbar = dbc.Navbar(
    dbc.Container([
        html.A(
            dbc.Row([
                dbc.Col(html.Img(src=app.get_asset_url("lapam.png"), className="lab-logo")),
                dbc.Col(dbc.NavbarBrand("Genomic Resistance Explorer", className="navbar-brand-text")),
            ], align="center", className="g-0"),
            href="/",
            style={"textDecoration": "none"},
        ),
        dbc.Nav([
            # Standard text link
            dbc.NavItem(dbc.NavLink("Home", href="/")),
            
            # Link with a custom image
            dbc.NavItem(
                dbc.NavLink(
                    [
                        html.Img(src=app.get_asset_url("code.png"), style={"height": "20px", "marginRight": "8px"}),
                    ], 
                    href="https://github.com/falatfernando/tbdashboard",
                    target="_blank",
                    className="d-flex align-items-center"
                )
            ),
            
            # Link with a custom image
            dbc.NavItem(
                dbc.NavLink(
                    [
                        html.Img(src=app.get_asset_url("github.png"), style={"height": "20px", "marginRight": "8px"}),
                    ], 
                    href="https://github.com/falatfernando", 
                    target="_blank",
                    className="d-flex align-items-center"
                )
            ),
        ], className="ms-auto", navbar=True)
    ], fluid=True),
    className="navbar-custom sticky-top mb-0"
)

# App layout
app.layout = html.Div([
    dcc.Store(id="current-gene-store"),
    navbar,
    
    # Hero / Search Section
    html.Div([
        dbc.Container([
            dbc.Row([
                dbc.Col([
                    html.H1([
                        dbc.Badge("TB", className="badge-tb me-2"), 
                        "Dashboard"
                    ], className="display-5 fw-bold mb-3"),
                    html.P(
                        "Explore M. tuberculosis WHO catalogue mutations and visualize the genomic region with precision.",
                        className="lead text-muted mb-5"
                    ),
                    
                    # Modern Search Card
                    dbc.Card([
                        dbc.CardBody([
                            dbc.InputGroup([
                                dbc.Input(
                                    id="search-input",
                                    type="text",
                                    placeholder="Search gene (e.g., rpoB) or mutation (e.g., katG_S315T)...",
                                    className="search-input",
                                    debounce=True
                                ),
                                dbc.Button(
                                    html.I(className="bi bi-search me-2"),
                                    id="search-button",
                                    className="search-btn"
                                )
                            ]),
                            html.Div([
                                html.Span("Quick access: ", className="text-muted small me-2"),
                                dbc.Button("rpoB", id="quick-rpoB", size="sm", className="quick-search-btn"),
                                dbc.Button("gyrA", id="quick-gyrA", size="sm", className="quick-search-btn"),
                                dbc.Button("katG", id="quick-katG", size="sm", className="quick-search-btn"),
                            ], className="mt-3 d-flex align-items-center justify-content-center")
                        ])
                    ], className="search-card")
                ], lg=8, className="mx-auto")
            ])
        ], className="py-5")
    ], className="hero-section"),
    
    # Results Section
    dbc.Container([
        dcc.Loading(
            id="loading-results",
            type="default",
            children=html.Div(id="search-results-container")
        ),
        
        # Footer
        html.Footer([
            dbc.Row([
                dbc.Col([
                    html.Hr(className="mb-4"),
                    html.Div([
                        html.Img(src=app.get_asset_url("lapam.png"), style={"height": "30px", "marginRight": "10px"}),
                        html.Span("TB Dashboard v1.1 | Developed by Fernando Falat", className="fw-bold")
                    ], className="d-flex align-items-center justify-content-center mb-3"),
                    html.P([
                        html.Strong("Data Sources: "),
                        "H37Rv (NC_000962.3), WHO Catalogue of mutations in Mycobacterium tuberculosis complex and their association with drug resistance, 2nd ed. | Visualization powered by JBrowse 2."
                    ], className="text-center small mb-0 text-muted")
                ])
            ])
        ], className="footer")
    ], fluid=False)
    
])


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
            return dbc.Alert([
                html.I(className="bi bi-exclamation-triangle-fill me-2"),
                f"Gene '{search_value}' not found. Please try another search."
            ], color="warning", className="rounded-3 shadow-sm"), None
    
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
                        html.H4(f"{gene_info.gene_name}", className="mb-0 fw-bold text-primary"),
                        html.Small(f"Locus: {gene_info.locus_tag}", className="text-muted")
                    ]),
                    dbc.Col([
                        html.Div([
                            html.Span("Product: ", className="fw-bold"),
                            html.Span(gene_info.product or 'Unknown')
                        ], className="text-end text-muted small")
                    ])
                ], align="center")
            ], className="card-header-custom"),
            dbc.CardBody([
                dbc.Row([
                    dbc.Col([
                        html.Div([
                            html.Label("Genomic Region", className="text-muted small d-block"),
                            html.Span(f"{gene_info.chromosome}:{gene_info.start:,}-{gene_info.end:,}", className="fw-500")
                        ])
                    ], width=4),
                    dbc.Col([
                        html.Div([
                            html.Label("Strand", className="text-muted small d-block"),
                            dbc.Badge("+" if gene_info.strand == "+" else "-", className="badge-custom")
                        ])
                    ], width=4),
                    dbc.Col([
                        html.Div([
                            html.Label("Length", className="text-muted small d-block"),
                            html.Span(f"{gene_length:,} bp", className="fw-500")
                        ])
                    ], width=4)
                ])
            ])
        ], className="result-card")
    )
    
    # JBrowse Card
    jbrowse_start = max(1, gene_info.start - 500)
    jbrowse_end = gene_info.end + 500
    jb_config = data_loader.get_jbrowse_config(gene_info.chromosome, jbrowse_start, jbrowse_end)
    
    results.append(
        dbc.Card([
            dbc.CardHeader([
                html.H5([html.I(className="bi bi-eye me-2"), "Genomic Visualization"], className="card-header-title")
            ], className="card-header-custom"),
            dbc.CardBody([
                html.Div([
                    dash_jbrowse.LinearGenomeView(
                        id="jbrowse-linear-view",
                        assembly=jb_config["assembly"],
                        tracks=jb_config["tracks"],
                        defaultSession=jb_config["defaultSession"],
                        location=f"{gene_info.chromosome}:{jbrowse_start}-{jbrowse_end}"
                    )
                ], className="jbrowse-container"),
                dbc.Button(
                    [html.I(className="bi bi-box-arrow-up-right me-2"), "Explore Full Genome"],
                    id="open-jbrowse",
                    color="link",
                    className="mt-3 p-0 text-decoration-none small",
                    href=f"#jbrowse?region={gene_info.chromosome}:{jbrowse_start}-{jbrowse_end}"
                )
            ])
        ], className="result-card")
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
                        dbc.Row([
                            dbc.Col(html.H6(drug, className="mb-0 fw-bold"), width="auto"),
                            dbc.Col(dbc.Badge(f"{len(drug_data)} mutations", className="badge-custom ms-2"), width="auto")
                        ], align="center")
                    ], className="bg-light border-0 py-2"),
                    dbc.CardBody([
                        html.P([
                            html.Strong("Evidence Tiers: ", className="small text-muted"),
                            html.Span(", ".join([str(t) for t in tiers if pd.notna(t)]), className="small fw-bold")
                        ], className="mb-3"),
                        html.Div(
                            id=f"drug-detail-{drug}",
                            children=create_drug_detail_table(drug_data, gene_info.gene_name)
                        )
                    ], className="pt-0")
                ], className="mb-3 border-0 shadow-sm")
            )
        
        results.append(
            dbc.Card([
                dbc.CardHeader([
                    html.H5([html.I(className="bi bi-capsule me-2"), "Drug Resistance Profile"], className="card-header-title")
                ], className="card-header-custom"),
                dbc.CardBody(drug_cards)
            ], className="result-card")
        )
    
    # Genomic Coordinates Card
    if len(mutations_df) > 0:
        results.append(
            dbc.Card([
                dbc.CardHeader([
                    dbc.Row([
                        dbc.Col(html.H5([html.I(className="bi bi-map me-2"), "Genomic Coordinates"], className="card-header-title")),
                        dbc.Col(html.Small(f"{len(mutations_df)} variants found", className="text-muted text-end"), width="auto")
                    ], align="center")
                ], className="card-header-custom"),
                dbc.CardBody([
                    html.P("Select a variant row below to view detailed coordinate calculations.", className="text-muted small mb-3"),
                    create_genomic_coords_table(mutations_df, gene_info)
                ])
            ], className="result-card")
        )
    
    # Calculator Card
    results.append(
        dbc.Card([
            dbc.CardHeader([
                html.H5([html.I(className="bi bi-calculator me-2"), "Coordinate Analysis"], className="card-header-title")
            ], className="card-header-custom"),
            dbc.CardBody([
                html.Div(id="coordinate-calculator-display", children=[
                    html.Div([
                        html.I(className="bi bi-info-circle me-2"),
                        "Selection required: Click a row in the Genomic Coordinates table above."
                    ], className="text-muted text-center py-4")
                ])
            ])
        ], id="calculator-card", className="result-card mb-5")
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
        return html.Div("Select a mutation from the table above to see detailed coordinate calculations.", className="text-muted")
    
    # Get the row data from the current viewport (handles pagination/filtering)
    row_idx = active_cell['row']
    
    if row_idx >= len(viewport_data):
        return html.Div("Error: Selected row not found.", className="text-danger")
        
    row = viewport_data[row_idx]
    variant = row['variant']
    gene_name = gene_data['gene_name']
    
    # Use the SPECIFIC genomic position from the clicked row
    actual_genomic_pos = int(row['position']) if str(row['position']).isdigit() else None
    ref = row.get('reference_nucleotide', 'N/A')
    alt = row.get('alternative_nucleotide', 'N/A')
    
    from data_utils import GeneInfo
    rel_pos = coord_calculator.calculate_relative_position(actual_genomic_pos, 
                                                         GeneInfo(
                                                             gene_id="", gene_name=gene_name, locus_tag=gene_data['locus_tag'],
                                                             start=gene_data['start'], end=gene_data['end'], strand=gene_data['strand'],
                                                             product=gene_data['product'], chromosome=gene_data['chromosome']
                                                         )) if actual_genomic_pos else None
    
    # Get ALL drug resistance entries for this variant
    resistance_df = data_loader.get_drug_resistance_info(gene_name, variant)
    
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
                        html.Strong(res_row['drug'], className="text-dark"),
                        html.Br(),
                        html.Small(f"Tier {res_row.get('tier', 'N/A')}", className="text-muted")
                    ]), width=4),
                    dbc.Col(dbc.Badge(
                        res_row.get('FINAL CONFIDENCE GRADING', 'N/A'), 
                        className="badge-tb" if "Assoc w R" in str(res_row.get('FINAL CONFIDENCE GRADING', '')) else "badge-custom"
                    ), width=8, className="text-end"),
                ], className="py-2 border-bottom align-items-center g-0")
            ]))
        res_display = html.Div(res_items, className="mt-2")
    else:
        res_display = html.P("No drug resistance data found in catalogue for this variant.", className="text-muted small")

    # Create display
    return html.Div([
        dbc.Row([
            dbc.Col([
                html.H5(variant, className="text-primary fw-bold mb-3"),
            ])
        ]),
        dbc.Row([
            dbc.Col([
                html.Div([
                    html.Label("Genomic Position", className="text-muted small d-block"),
                    html.Span(gen_pos_display, className="fw-bold")
                ], className="mb-3"),
                html.Div([
                    html.Label("Nucleotide Change", className="text-muted small d-block"),
                    html.Span(f"{ref} > {alt}", className="fw-bold")
                ])
            ], width=6),
            dbc.Col([
                html.Div([
                    html.Label("Gene Relative (c.)", className="text-muted small d-block"),
                    html.Span(f"c.{rel_pos}" if rel_pos else "N/A", className="fw-bold text-success")
                ], className="mb-3"),
                html.Div([
                    html.Label("Gene Context", className="text-muted small d-block"),
                    html.Span(f"{gene_name} ({gene_data['strand']})", className="fw-bold")
                ])
            ], width=6)
        ], className="mb-4 bg-light p-3 rounded-3"),
        
        html.H6("Resistance Association", className="fw-bold mb-2 small text-uppercase letter-spacing-1"),
        res_display,
        
        dbc.Accordion([
            dbc.AccordionItem([
                html.Code(
                    f"Formula ({gene_data['strand']} strand):\n"
                    f"{'genomic - start + 1' if gene_data['strand'] == '+' else 'end - genomic + 1'}\n\n"
                    f"Calculation:\n"
                    f"{f'{actual_genomic_pos:,} - {gene_data['start']:,} + 1' if gene_data['strand'] == '+' else f'{gene_data['end']:,} - {actual_genomic_pos:,} + 1'} = {rel_pos}",
                    style={"whiteSpace": "pre-wrap", "display": "block", "padding": "15px", "fontSize": "0.85rem"}
                )
            ], title="View Mathematical Derivation", className="mt-4 border-0")
        ], start_collapsed=True, className="mt-4")
    ])


def create_drug_detail_table(drug_data: pd.DataFrame, gene_name: str) -> html.Div:
    """Create detailed drug resistance table."""
    tier_order = {'1': 1, '2': 2, '3': 3, '4': 4, '5': 5}
    drug_data_sorted = drug_data.copy()
    drug_data_sorted['tier_sort'] = drug_data_sorted['tier'].map(tier_order).fillna(99)
    drug_data_sorted = drug_data_sorted.sort_values('tier_sort')
    
    cols = [
        {'name': 'Mutation', 'id': 'mutation'},
        {'name': 'Tier', 'id': 'tier'},
        {'name': 'Confidence', 'id': 'FINAL CONFIDENCE GRADING'}
    ]
    
    table = dash_table.DataTable(
        data=drug_data_sorted.to_dict('records'),
        columns=cols,
        filter_action="native",
        sort_action="native",
        page_size=5,
        style_table={'overflowX': 'auto'},
        style_cell={
            'textAlign': 'left',
            'padding': '12px',
            'fontFamily': 'inherit',
            'fontSize': '0.9rem'
        },
        style_header={
            'backgroundColor': 'transparent',
            'fontWeight': 'bold',
            'borderBottom': '2px solid #E2E8F0',
            'textTransform': 'uppercase',
            'fontSize': '0.75rem',
            'letterSpacing': '0.05em'
        },
        style_data={
            'borderBottom': '1px solid #F1F5F9'
        }
    )
    
    return html.Div(table, className="dash-table-container")


def create_genomic_coords_table(mutations_df: pd.DataFrame, gene_info) -> html.Div:
    """Create genomic coordinates drill-down table."""
    display_df = mutations_df.copy()
    
    display_df['Gene Relative Pos'] = display_df.apply(
        lambda row: coord_calculator.calculate_relative_position(
            int(row['position']), gene_info
        ) if row.get('position') and str(row['position']).isdigit() else None,
        axis=1
    )
    
    display_df['Gene Relative'] = display_df['Gene Relative Pos'].apply(
        lambda x: f"c.{int(x)}" if pd.notna(x) else 'N/A'
    )
    
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
        cell_selectable=True,
        filter_action="native",
        sort_action="native",
        page_size=10,
        style_table={'overflowX': 'auto'},
        style_cell={
            'textAlign': 'left',
            'padding': '12px',
            'fontFamily': 'inherit',
            'fontSize': '0.9rem'
        },
        style_header={
            'backgroundColor': '#F8FAFC',
            'fontWeight': 'bold',
            'borderBottom': '2px solid #E2E8F0'
        },
        style_data_conditional=[
            {
                'if': {'state': 'active'},
                'backgroundColor': 'rgba(156, 203, 203, 0.1)',
                'border': '1px solid #438E8E'
            }
        ]
    )
    
    return html.Div(table, className="dash-table-container")


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=8050)
