import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pickle
import numpy as np
import ribopy
from ribopy import Ribo
from functions import get_cds_range_lookup, get_sequence
from functions_filter import get_filtered_transcripts, get_filtered_zscores
from functions_heatmap_v4 import find_common_stall_sites, collect_stall_sequences, create_raw_heatmap, normalize_heatmap, get_heatmap_df, get_stall_sites_df
import base64
import io
import tempfile
import traceback
import os
import pandas as pd

app = dash.Dash(__name__)

app.layout = html.Div([
    html.H1("Heatmap", style={'margin-left': '25px'}),
    html.Div([
        html.Div([
            html.Label('Upload ribo file:', style={'margin-left': '25px'}),
            dcc.Upload(
                id='upload-ribo-file',
                children=html.Div(['Drag and Drop or ', html.A('Select File')]),
                style={
                    'width': '80%',
                    'height': '60px',
                    'lineHeight': '60px',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'margin-left': '20px'
                },
                multiple=False
            )
        ], style={'flex': '1', 'margin-right': '20px'}),
        html.Div([
            html.Label('Upload pickle file:', 
                        style={'margin-left': '25px'}),
            dcc.Upload(
                id='upload-pickle-file',
                children=html.Div(['Drag and Drop or ', html.A('Select File')]),
                style={
                    'width': '80%',
                    'height': '60px',
                    'lineHeight': '60px',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'margin-left': '20px'
                },
                multiple=False
            )
        ], style={'flex': '1', 'margin-right': '20px'}),
        html.Div([
            html.Label('Upload reference file:', 
                        style={'margin-left': '25px'}),
            dcc.Upload(
                id='upload-reference-file',
                children=html.Div(['Drag and Drop or ', html.A('Select File')]),
                style={
                    'width': '80%',
                    'height': '60px',
                    'lineHeight': '60px',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'margin-left': '20px'
                },
                multiple=False
            )
        ], style={'flex': '1'})
    ], style={'display': 'flex', 'justify-content': 'space-between', 'margin-left': '25px'}),
    html.Div([
        html.Label('Select organism: ', 
                    style={'margin-left': '25px'}),
        dcc.Dropdown(
            id='organism',
            options=[
                {'label': 'Mouse', 'value': 1},
                {'label': 'Other', 'value': 2},
            ],
            value=1,
            style={'width': '99%', 'margin-left': '10px'}
        ),
    ], style={'margin-top': '20px'}),
    html.Div([
        html.Label('Enter a list of replicates of an experimental condition, separated by a space. To look at multiple conditions, separate by new line.',
                    style={'margin-left': '25px'}),
        html.P('Example input:', 
                style={'margin-left': '50px', 'padding': '0'}),
        html.P('Control_1 Control_2 Control_3',
                style={'margin-left': '50px', 'line-height': '0.25'}),
        html.P('Mutant_1 Mutant_2 Mutant_3', 
                style={'margin-left': '50px', 'line-height': '0.1'}),
        dcc.Textarea(id='experiments', style={'margin-left': '25px', 'width': '97%', 'height': '100px'}),
    ], style={'margin-top': '20px'}),
    html.Div([
        html.Label('Enter number of transcripts to use, determined by the number of reads per nucleotide within the coding region: ', 
                style={'margin-left': '25px'}),
        dcc.Input(id='num_transcripts', type='number', value=100, style={'width': '4%'}),
        html.P('If there is more than one experiment, the intersection of transcripts will be used.',
            style={'margin-left': '25px', 'margin-top': '0px'}),
    ], style={'margin-top': '20px'}), 
    html.Div([
        html.Label('Enter percentile threshold for determining stall sites (0 to 100): ',
                    style={'margin-left': '25px'}),
        dcc.Input(id='percentile', type='number', value=99, style={'width': '3%'}),
    ], style={'margin-top': '20px'}),
    html.Button('Generate Heatmap', id='generate_button', n_clicks=0, style={'margin-left': '25px', 'margin-top': '20px'}),
    dcc.Loading(
        id="loading",
        children=[html.Div(id='heatmap_output')],
        type="default",
    )
])

@app.callback(
    Output('heatmap_output', 'children'),
    Input('generate_button', 'n_clicks'),
    State('upload-ribo-file', 'contents'),
    State('upload-ribo-file', 'filename'),
    State('upload-pickle-file', 'contents'),
    State('upload-pickle-file', 'filename'),
    State('upload-reference-file', 'contents'),
    State('upload-reference-file', 'filename'),
    State('organism', 'value'),
    State('experiments', 'value'),
    State('num_transcripts', 'value'),
    State('percentile', 'value')
)
def update_heatmap(n_clicks, ribo_content, ribo_filename, pickle_content, pickle_filename, ref_content, ref_filename, organism, experiments_str, num_transcripts, percentile):
    if n_clicks > 0:
        # Set alias based on organism
        alias = True if organism == 1 else False
        
        # Parse experiments input
        experiments = []
        for replicates_str in experiments_str.split('\n'):
            if replicates_str:
                replicates = replicates_str.split(' ')
                experiments.append(replicates)
        
        # Decode ribo file
        content_type, content_string = ribo_content.split(',')
        decoded = base64.b64decode(content_string)
        try:
            # Create a temporary file for ribo file
            with tempfile.NamedTemporaryFile(delete=False) as temp_file:
                temp_file.write(decoded)
                temp_file_path = temp_file.name
            
            if alias == True:
                ribo_object = Ribo(temp_file_path, alias=ribopy.api.alias.apris_human_alias)
            else: 
                ribo_object = Ribo(temp_file_path)
            cds_range = get_cds_range_lookup(ribo_object)
        except Exception as e:
            traceback.print_exc()
            return html.Div('Invalid ribo file.', style={'color': 'red'})

        try:
            pickle_data = base64.b64decode(pickle_content.split(',')[1])
            with tempfile.NamedTemporaryFile(delete=False, suffix=os.path.splitext(pickle_filename)[1]) as pickle_temp_file:
                pickle_temp_file.write(pickle_data)
                pickle_file_path = pickle_temp_file.name
        except Exception as e:
            traceback.print_exc()
            return html.Div('Invalid pickle file.', style={'color': 'red'})
        try:
            # Decode reference file and create a temporary file
            ref_data = base64.b64decode(ref_content.split(',')[1])
            with tempfile.NamedTemporaryFile(delete=False, suffix=os.path.splitext(ref_filename)[1]) as ref_temp_file:
                ref_temp_file.write(ref_data)
                reference_file_path = ref_temp_file.name

            # Get sequence and CDS range lookup
            sequence = get_sequence(ribo_object, reference_file_path, alias)
            cds_range = get_cds_range_lookup(ribo_object)
            transcripts = get_filtered_transcripts(pickle_file_path, list(np.array(experiments).flat), num_transcripts)

            # Main loop to generate heatmaps
            output_excel_file = 'codon_heatmaps.xlsx'
            with pd.ExcelWriter(output_excel_file, engine='xlsxwriter') as writer:
                all_norm_heatmaps = []
                for replicates in experiments:
                    common_stall_sites = find_common_stall_sites(replicates, pickle_file_path, transcripts, percentile)
                    stall_sequences = collect_stall_sequences(common_stall_sites, sequence, cds_range)
                    raw_heatmap = create_raw_heatmap(stall_sequences)
                    norm_heatmap = normalize_heatmap(raw_heatmap, sequence, cds_range, transcripts)
                    all_norm_heatmaps.append(norm_heatmap)

                    df_heatmap = get_heatmap_df(raw_heatmap, sequence, cds_range, transcripts)
                    df_heatmap.to_excel(writer, sheet_name=f'heatmap_{replicates[0]}')
                    df_stallsites = get_stall_sites_df(common_stall_sites, cds_range, sequence)
                    df_stallsites.to_excel(writer, sheet_name=f'stall_sites_{replicates[0]}')

            # Calculate global zmin and zmax
            all_z_values = np.concatenate([norm_heatmap.values.flatten() for norm_heatmap in all_norm_heatmaps])
            zmin, zmax = np.min(all_z_values), np.max(all_z_values)

            # Create a subplot figure
            fig = make_subplots(
                rows=1, cols=len(all_norm_heatmaps),
                subplot_titles=[f"Heatmap for {replicates[0]}" for replicates in experiments]
            )

            # Add each heatmap to the subplot
            for i, norm_heatmap in enumerate(all_norm_heatmaps):
                heatmap = go.Heatmap(
                    z=norm_heatmap.values,
                    x=norm_heatmap.columns,
                    y=norm_heatmap.index,
                    colorscale='Viridis',
                    zmin=zmin,
                    zmax=zmax
                )
                fig.add_trace(heatmap, row=1, col=i+1)

            # Update layout settings for each subplot
            for i in range(len(all_norm_heatmaps)):
                fig.update_xaxes(tickfont=dict(size=10), dtick=1, row=1, col=i+1)
                fig.update_yaxes(tickfont=dict(size=8), row=1, col=i+1)

            fig.update_layout(
                height=850,
                width=(1400/len(all_norm_heatmaps)) * len(all_norm_heatmaps),
                title_text="Combined Heatmaps"
            )

            # Return the figure to the output div
            return dcc.Graph(figure=fig)

        except Exception as e:
            # Log the full traceback for debugging
            traceback.print_exc()
            return html.Div('Error processing reference file.', style={'color': 'red'})
        finally:
            # Ensure the temporary files are deleted
            try:
                os.remove(temp_file_path)
            except Exception as e:
                print(f"Error deleting temporary ribo file: {e}")
            try:
                os.remove(pickle_file_path)
            except Exception as e:
                print(f"Error deleting temporary pickle file: {e}")
            try:
                os.remove(reference_file_path)
            except Exception as e:
                print(f"Error deleting temporary reference file: {e}")

    return html.Div()

if __name__ == '__main__':
    app.run_server(debug=True)