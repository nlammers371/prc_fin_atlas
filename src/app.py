import dash
import plotly.express as px
from dash import dcc
from dash import html
from dash.dependencies import Input, Output, State
import dash_daq as daq
import plotly.graph_objs as go
import pandas as pd
import numpy as np
import json
import glob2 as glob
import os
from scipy.interpolate import LinearNDInterpolator
import alphashape

########################
# App
app = dash.Dash(__name__)  # , external_stylesheets=external_stylesheets)
server = app.server

# global init_toggle
# init_toggle = True

dataRoot = "../dat/"


global fileList, imNameList
# get list of filepaths
fileList = sorted(glob.glob(dataRoot + '*.csv'))
imNameList = []
for fn in range(len(fileList)):
    labelName = fileList[fn].replace(dataRoot, '', 1)
    labelName = labelName.replace('.csv', '')
    imNameList.append(labelName)

# compile dictionary of gene names that correspond to each dataset
gene_name_dict = {}
for f in range(len(fileList)):
    filename_temp = imNameList[f]
    # propPath = dataRoot + filename_temp + '_nucleus_props.csv'

    # if os.path.isfile(propPath):
    # df_temp = pd.read_csv(propPath, index_col=0)
    # else:
    #     raise Exception(
    #         f"Selected dataset( {filename_temp} ) dataset has no nucleus data. Have you run extract_nucleus_stats?")
    #
    # # get list of genes that we can look at
    # colnames_temp = df_temp.columns
    # list_raw = [item for item in colnames_temp if "_cell_mean_nn" in item]
    # gene_names_temp = [item.replace("_cell_mean_nn", "") for item in list_raw]

    gene_name_dict[filename_temp] = ['Col11a2', 'Emilin3a', 'Fgf10a', 'Hand2', 'Myod1', 'Prdm1a', 'Robo3', 'Sox9a',
                                     'Tbx5a']

def load_nucleus_dataset(filename):
    # global fin_points_prev, not_fin_points_prev, class_predictions_curr, df, curationPath, propPath

    propPath = dataRoot + filename + '.csv'

    # if os.path.isfile(propPath):
    df_raw = pd.read_csv(propPath, index_col=0)

    if "_hr" in filename:
        keep_indices = np.where((df_raw["Z"] <= 0.46))
        df = df_raw.iloc[keep_indices].copy()
    else:
        df = df_raw
    # else:
    #     raise Exception(
    #         f"Selected dataset( {filename} ) dataset has no nucleus data. Have you run extract_nucleus_stats?")

    # fin_nuclei = np.where(df["pec_fin_flag"] == 2)
    # df = df.iloc[fin_nuclei]

    # normalize gene expression levels
    # colnames = df.columns
    # list_raw = [item for item in colnames if "_cell_mean_nn" in item]
    # gene_names = [item.replace("_cell_mean_nn", "") for item in list_raw]

    gene_names = ['Col11a2', 'Emilin3a', 'Fgf10a', 'Hand2', 'Myod1', 'Prdm1a', 'Robo3', 'Sox9a', 'Tbx5a']

    for g in range(len(gene_names)):
        # ind_list = [i for i in range(len(colnames)) if g in colnames[i]]

        colname = gene_names[g]
        c_max = np.max(df[colname])
        df[colname] = df.loc[:, colname] / c_max

    return df

global df, colnames, gene_names

df = load_nucleus_dataset(imNameList[0])

# get list of genes that we can look at
# colnames = df.columns
# list_raw = [item for item in colnames if "_cell_mean_nn" in item]
gene_names = ['Col11a2', 'Emilin3a', 'Fgf10a', 'Hand2', 'Myod1', 'Prdm1a', 'Robo3', 'Sox9a', 'Tbx5a']

plot_list = ["3D Scatter", "Volume Plot"]


def create_figure(df, gene_name=None, plot_type=None):
    colormaps = px.colors.named_colorscales()  # ["ice", "inferno", "viridis"]
    colormaps[2] = "inferno"
    colormaps[3] = "matter"

    if gene_name == None:
        plot_gene = gene_names[2]  # + "_mean_nn"
        cmap = colormaps[2]
    else:
        plot_gene = gene_name  # + "_mean_nn"
        g_index = gene_names.index(gene_name)
        cmap = colormaps[g_index]

    xs = "X"
    ys = "Y"
    zs = "Z"
    df = df.dropna(subset=["Z"])
    alph = 4

    xyz_array = np.asarray(df[[xs, ys, zs]])

    if (plot_type == None) | (plot_type == "3D Scatter"):
        high_flags = df[plot_gene] >= 0.3
        low_flags = df[plot_gene] < 0.3
        fig = px.scatter_3d(df.iloc[np.where(high_flags)], x=xs, y=ys, z=zs, opacity=0.8, color=plot_gene,
                            color_continuous_scale=cmap, range_color=(0, 1))

        fig.update_traces(marker=dict(size=8
                                      ))

        fig.add_trace(go.Scatter3d(x=df[xs].iloc[np.where(low_flags)],
                                   y=df[ys].iloc[np.where(low_flags)],
                                   z=df[zs].iloc[np.where(low_flags)],
                                   mode='markers',
                                   marker=dict(
                                       size=5,
                                       color=df[plot_gene].iloc[np.where(low_flags)],
                                       # set color to an array/list of desired values
                                       colorscale=cmap,  # choose a colorscale
                                       opacity=0.3,
                                       cmin=0,
                                       cmax=1)
                                   ))
        fig.update_layout(coloraxis_colorbar_title_text='Normalized Expression')
        fig.update_layout(showlegend=False)
        # fig = px.scatter_3d(df.iloc[high_flags], x=xs, y=ys, z=zs, opacity=0.6, size=plot_gene, color=plot_gene, color_continuous_scale=cmap)
        # fig.update_traces(marker=dict(size=8))

        # fig.add_trace(go.Mesh3d(x=xyz_array[:, 0], y=xyz_array[:, 1], z=xyz_array[:, 2],
        #                         alphahull=alph,
        #                         opacity=0.1,
        #                         color='gray'))

    elif (plot_type == "Volume Plot"):

        # generate points to interpolate
        xx = np.linspace(min(df[xs]), max(df[xs]), num=30)
        yy = np.linspace(min(df[ys]), max(df[ys]), num=30)
        zz = np.linspace(min(df[zs]), max(df[zs]), num=30)

        X, Y, Z = np.meshgrid(xx, yy, zz)  # 3D grid for interpolation

        # generate normalized arrays
        X_norm = X / np.max(X)
        Y_norm = Y / np.max(Y)
        Z_norm = Z / np.max(Z)

        # generate interpolator
        gene_values = np.asarray(df[plot_gene])
        interp = LinearNDInterpolator(xyz_array, gene_values)

        # get interpolated estimate of gene expression
        G = interp(X, Y, Z)

        # generate alpha shape
        xyz_array_norm = np.divide(xyz_array, np.asarray([np.max(X), np.max(Y), np.max(Z)]))
        alpha_fin = alphashape.alphashape(xyz_array_norm, 9)
        xyz_long = np.concatenate((np.reshape(X_norm, (X.size, 1)),
                                   np.reshape(Y_norm, (Y.size, 1)),
                                   np.reshape(Z_norm, (Z.size, 1))),
                                  axis=1)
        inside_flags = alpha_fin.contains(xyz_long)
        G_long = G.flatten()
        G_long[~inside_flags] = np.nan

        fig = go.Figure(
            data=go.Volume(
                x=X.flatten(),
                y=Y.flatten(),
                z=Z.flatten(),
                value=G_long,
                opacity=.1,
                isomin=0.2,
                isomax=0.8,  # needs to be small to see through all surfaces
                surface_count=15,  # needs to be a large number for good volume rendering
                colorscale=cmap
            ))

        fig.add_trace(go.Mesh3d(x=xyz_array[:, 0], y=xyz_array[:, 1], z=xyz_array[:, 2],
                                alphahull=4,
                                opacity=0.1,
                                color='gray'))

        fig.update_layout(coloraxis_colorbar_title_text='Normalized Expression')

    elif plot_type == "Multiplot":
        # fig = px.scatter_3d(df, x=xs, y=ys, z=zs, opacity=0.4)
        # raise Warning("Plot type " + plot_type + " is not currently supported")
        # generate points to interpolate
        xx = np.linspace(min(df[xs]), max(df[xs]), num=25)
        yy = np.linspace(min(df[ys]), max(df[ys]), num=25)
        zz = np.linspace(min(df[zs]), max(df[zs]), num=25)

        X, Y, Z = np.meshgrid(xx, yy, zz)  # 3D grid for interpolation

        # generate normalized arrays
        X_norm = X / np.max(X)
        Y_norm = Y / np.max(Y)
        Z_norm = Z / np.max(Z)

        # generate interpolator
        gene_values1 = np.asarray(df[gene_names[0] + "_mean_nn"])
        interp1 = LinearNDInterpolator(xyz_array, gene_values1)

        gene_values2 = np.asarray(df[gene_names[1] + "_mean_nn"])
        interp2 = LinearNDInterpolator(xyz_array, gene_values2)

        gene_values3 = np.asarray(df[gene_names[2] + "_mean_nn"])
        interp3 = LinearNDInterpolator(xyz_array, gene_values3)

        # get interpolated estimate of gene expression
        G1 = interp1(X, Y, Z)
        G2 = interp2(X, Y, Z)
        G3 = interp3(X, Y, Z)

        # generate alpha shape
        xyz_array_norm = np.divide(xyz_array, np.asarray([np.max(X), np.max(Y), np.max(Z)]))
        alpha_fin = alphashape.alphashape(xyz_array_norm, 9)
        xyz_long = np.concatenate((np.reshape(X_norm, (X.size, 1)),
                                   np.reshape(Y_norm, (Y.size, 1)),
                                   np.reshape(Z_norm, (Z.size, 1))),
                                  axis=1)

        xyz_long2 = np.concatenate((np.reshape(X, (X.size, 1)),
                                    np.reshape(Y, (Y.size, 1)),
                                    np.reshape(Z, (Z.size, 1))),
                                   axis=1)

        inside_flags = alpha_fin.contains(xyz_long)

        G1_long = G1.flatten()
        G1_long[~inside_flags] = np.nan
        G1_indices = np.where(G1_long >= 0.4)

        G2_long = G2.flatten()
        G2_long[~inside_flags] = np.nan
        G2_indices = np.where(G2_long >= 0.4)

        G3_long = G3.flatten()
        G3_long[~inside_flags] = np.nan
        G3_indices = np.where(G3_long >= 0.4)

        fig = go.Figure(data=go.Mesh3d(x=xyz_long2[G1_indices, 0][0], y=xyz_long2[G1_indices, 1][0],
                                       z=xyz_long2[G1_indices, 2][0],
                                       alphahull=9,
                                       opacity=0.2,
                                       color='blue',
                                       name=gene_names[0]))

        # print(xyz_long2[G1_indices, 0])
        fig.add_trace(
            go.Mesh3d(x=xyz_long2[G2_indices, 0][0], y=xyz_long2[G2_indices, 1][0], z=xyz_long2[G2_indices, 2][0],
                      alphahull=15,
                      opacity=0.2,
                      color='red',
                      name=gene_names[1]))

        fig.add_trace(
            go.Mesh3d(x=xyz_long2[G3_indices, 0][0], y=xyz_long2[G3_indices, 1][0], z=xyz_long2[G3_indices, 2][0],
                      alphahull=9,
                      opacity=0.2,
                      color='green',
                      name=gene_names[2]))

        fig.add_trace(go.Mesh3d(x=xyz_array[:, 0], y=xyz_array[:, 1], z=xyz_array[:, 2],
                                alphahull=9,
                                opacity=0.1,
                                color='gray'))

        fig.update_layout(showlegend=True)
        fig.update_layout(legend_title_text='Gene Name')

    else:
        fig = px.scatter_3d(df, x=xs, y=ys, z=zs, opacity=0.4)
        raise Warning("Plot type " + plot_type + " is not currently supported")

    fig.update_layout(scene=dict(
        xaxis_title='A-P',
        yaxis_title='D-V',
        zaxis_title='P-D'))
    return fig

f = create_figure(df)

app.layout = html.Div([
    dcc.Graph(id='3d_scat', figure=f, style={'width': '90vh', 'height': '90vh', 'padding-left': '450px',
                                             'justify': 'center'}),

    html.Div(id='df_list', hidden=True),
    html.Div([
        html.Label(['Dataset:'], style={'font-weight': 'bold', "text-align": "center"}),
        dcc.Dropdown(options=[
                            {'label': 'Low-res atlas (250 cells)', 'value': imNameList[0]},
                            {'label': 'High-res atlas (1000 cells)', 'value': imNameList[1]},
                        ],
                     value=imNameList[1], id='dataset-dropdown'),
    ],
        style={'width': '33%', 'display': 'inline-block'}),
    html.Div(id='dd-output-container', hidden=True),

    html.Div([
        html.Label(['Gene:'], style={'font-weight': 'bold', "text-align": "center"}),
        dcc.Dropdown(id='gene-dropdown'),
    ],
        style={'width': '33%', 'display': 'inline-block'}),
    html.Div(id='gene-output-container', hidden=True),

    # html.Div([
    #     dcc.Dropdown(plot_list, id='plot-dropdown'),
    #     html.Div(id='plot-output-container')
    # ])
    html.Div(id='plot_list', hidden=True),
    html.Div([
        html.Label(['Plot Type:'], style={'font-weight': 'bold', "text-align": "center"}),
        dcc.Dropdown(plot_list, plot_list[0], id='plot-dropdown'),
    ],
        style={'width': '33%', 'display': 'inline-block'}),
    html.Div(id='plot-output-container', hidden=True),
    # html.Div(
    #     daq.ToggleSwitch(
    #         id='coord-toggle-switch',
    #         value=False,
    #         color='lightgreen'
    #     )),
    # html.Div(id='my-toggle-switch-output')
])

# @app.callback(
#     Output('my-toggle-switch-output', 'children'),
#     Input('coord-toggle-switch', 'value'))
# def update_output(value):
#     value = json.dumps(value)
#     toggle_val = value == "true"
#     if toggle_val:
#         toggle_string = "native"
#     else:
#         toggle_string = "unfolded"
#     return f'Showing {toggle_string} nucleus coordinates.'

# update list of genes
@app.callback(
    [dash.dependencies.Output('gene-dropdown', 'options'),
     dash.dependencies.Output('gene-dropdown', 'value')],
    [dash.dependencies.Input('dataset-dropdown', 'value')]
)
def update_date_dropdown(name):
    new_dict = [{'label': i, 'value': i} for i in gene_name_dict[name]]
    return new_dict, new_dict[2]['value']

# @app.callback(
#     Output('dd-output-container', 'children'),
#     Input('dataset-dropdown', 'value')
# )
# def load_wrapper(value):
#     return value
#
# @app.callback(
#     Output('gene-output-container', 'children'),
#     Input('gene-dropdown', 'value')
# )
# def load_wrapper(value):
#     return value

@app.callback(Output('3d_scat', 'figure'),
              [Input('gene-dropdown', 'value'),
               Input('dataset-dropdown', 'value'),
               Input('plot-dropdown', 'value')
               # Input('my-toggle-switch-output', 'children')
               ])
def chart_3d(gene_name, fileName, plot_type):

    global f, df, gene_names

    # check to see which values have changed
    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
    # toggle_bool = toggle_switch == "Showing unfolded nucleus coordinates."
    df = load_nucleus_dataset(fileName)

    # get list of genes that we can look at
    # colnames = df.columns
    # list_raw = [item for item in colnames if "_cell_mean_nn" in item]
    gene_names = ['Col11a2', 'Emilin3a', 'Fgf10a', 'Hand2', 'Myod1', 'Prdm1a', 'Robo3', 'Sox9a',
                  'Tbx5a']  # [item.replace("_cell_mean_nn", "") for item in list_raw]

    f = create_figure(df, gene_name, plot_type)

    f.update_layout(uirevision="Don't change")
    return f



if __name__ == '__main__':
    # app.run_server(debug=True)
    # app.run_server(debug=True, server=app.server)
    # noinspection PyPackageRequirements
    app.run_server(debug=True, port=8056)
    # set parameters


    # # load image data
    # visualize_atlas(dataRoot)