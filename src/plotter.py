import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def create_combined_plots(df: pd.DataFrame, xScatter: str, yScatter: str) -> None:
    # Create subplots with 2 rows
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=("Borghi Diagram", "Design Viability Map"),
        vertical_spacing=0.12,
        row_heights=[0.5, 0.5]
    )
    
    # ==================== BORGHI DIAGRAM (Top) ====================
    
    # Laminar curve
    x_lam = np.array([0.1, 10])
    y_lam = np.array([10, 0.1])
    fig.add_trace(go.Scatter(
        x=x_lam, y=y_lam, mode="lines",
        line=dict(color="black", width=1),
        name="Re_t = 1",
        showlegend=False
    ), row=1, col=1)
    
    fig.add_annotation(
        x=np.log10(0.5), y=np.log10(0.3), text="Laminar",
        showarrow=False, font=dict(size=14),
        xref="x1", yref="y1"
    )
    
    # u_prime/S_L = 1
    x_2 = np.array([1.0, 1000])
    y_2 = np.array([1.0, 1])
    fig.add_trace(go.Scatter(
        x=x_2, y=y_2, mode="lines",
        line=dict(color="black", width=1),
        name="u'/S_L = 1",
        showlegend=False
    ), row=1, col=1)
    
    fig.add_annotation(
        x=np.log10(400), y=np.log10(0.85), text="u'/S_L = 1",
        showarrow=False, font=dict(size=14),
        xref="x1", yref="y1"
    )
    
    # Ka = 1
    x_3 = np.linspace(1.0, 1000, 1000)
    y_3 = x_3**(1/3)
    fig.add_trace(go.Scatter(
        x=x_3, y=y_3, mode="lines",
        line=dict(color="black", width=1),
        name="Ka = 1",
        showlegend=False
    ), row=1, col=1)
    
    fig.add_annotation(
        x=np.log10(400), y=np.log10(10), text="Ka = 1",
        showarrow=False, font=dict(size=14),
        xref="x1", yref="y1"
    )
    
    # Ka = 10
    x_3a = np.linspace(0.317, 1000, 1000)
    y_3a = (10 * (x_3a)**(1/2))**(2/3)
    fig.add_trace(go.Scatter(
        x=x_3a, y=y_3a, mode="lines",
        line=dict(color="black", width=1, dash="dash"),
        name="Ka = 10",
        showlegend=False
    ), row=1, col=1)
    
    # Ka = 100
    x_3b = np.linspace(0.1, 1000, 1000)
    y_3b = (100 * (x_3b)**(1/2))**(2/3)
    fig.add_trace(go.Scatter(
        x=x_3b, y=y_3b, mode="lines",
        line=dict(color="black", width=1),
        name="Ka = 100",
        showlegend=False
    ), row=1, col=1)
    
    fig.add_annotation(
        x=np.log10(2), y=np.log10(36), text="Ka = 100",
        showarrow=False, font=dict(size=14),
        xref="x1", yref="y1"
    )
    
    # Da = 1
    x_4 = np.linspace(1.0, 1000, 1000)
    y_4 = x_4
    fig.add_trace(go.Scatter(
        x=x_4, y=y_4, mode="lines",
        line=dict(color="black", width=1),
        name="Da = 1",
        showlegend=False
    ), row=1, col=1)
    
    fig.add_annotation(
        x=np.log10(90), y=np.log10(150), text="Da = 1",
        showarrow=False, font=dict(size=14),
        xref="x1", yref="y1"
    )
    
    # Scatter data for Borghi
    df['hover_text'] = df.apply(lambda row: f"U0: {row['U0']}<br>PHI: {row['PHI']}", axis=1)
    
    for stab_type in df['Stabilization'].unique():
        df_subset = df[df['Stabilization'] == stab_type]
        fig.add_trace(go.Scatter(
            x=df_subset[xScatter],
            y=df_subset[yScatter],
            mode='markers',
            marker=dict(color="red", size=16, symbol='star'),
            name=stab_type,
            hovertext=df_subset['hover_text'],
            legendgroup='borghi'
        ), row=1, col=1)
    
    # ==================== DESIGN MAP (Bottom) ====================
    
    for stab_type in df['Stabilization'].unique():
        df_subset = df[df['Stabilization'] == stab_type]
        fig.add_trace(go.Scatter(
            x=df_subset['dP_percent'],
            y=df_subset['Da'],
            mode='markers',
            marker=dict(color="red", size=16, symbol='star'),
            name=stab_type,
            hovertemplate='<b>U0</b>: %{customdata[0]}<br>' +
                          '<b>PHI</b>: %{customdata[1]}<br>' +
                          '<b>L_premixer</b>: %{customdata[2]}<br>' +
                          '<b>BR</b>: %{customdata[3]}<br>' +
                          '<b>dP</b>: %{x}%<br>' +
                          '<b>Da</b>: %{y}<extra></extra>',
            customdata=df_subset[['U0', 'PHI', 'L_premixer', 'BR']].values,
            legendgroup='design',
            showlegend=False
        ), row=2, col=1)
    
    # Add constraint lines to design map
    fig.add_vline(
        x=2.0, line_width=3, line_dash="dash", line_color="red",
        row=2, col=1
    )
    fig.add_annotation(
        x=2.0, y=df['Da'].max()*0.95,
        text="Max Pressure Loss (2%)",
        showarrow=False, font=dict(size=12),
        xref="x2", yref="y2", xanchor="right"
    )
    
    fig.add_hline(
        y=1.0, line_width=3, line_dash="dash", line_color="black",
        row=2, col=1
    )
    fig.add_annotation(
        x=3.8, y=1.0,
        text="Da = 1",
        showarrow=False, font=dict(size=12),
        xref="x2", yref="y2", yanchor="bottom"
    )
    
    # Safe zone rectangle
    fig.add_shape(
        type="rect",
        x0=0, y0=1.0, x1=2.0, y1=df['Da'].max()*1.05,
        line=dict(color="Green", width=2),
        fillcolor="Green", opacity=0.1,
        row=2, col=1
    )
    
    # Update axes
    fig.update_xaxes(
        type="log",
        range=[-1, 3],
        showgrid=False, showline=True, mirror=True,
        title=r"$l_t/\delta_L$",
        row=1, col=1
    )
    fig.update_yaxes(
        type="log",
        range=[-1, 2.8],
        showgrid=False, showline=True, mirror=True,
        title=r"$u'/S_L$",
        row=1, col=1
    )
    
    fig.update_xaxes(
        range=[0, 4],
        showgrid=False, showline=True, mirror=True,
        title="Pressure Drop (%)",
        row=2, col=1
    )
    fig.update_yaxes(
        showgrid=False, showline=True, mirror=True,
        title="Damköhler Number (Da)",
        row=2, col=1
    )
    
    # Update layout
    fig.update_layout(
        template="simple_white",
        width=600,
        height=800,
        margin=dict(l=20, r=20, t=60, b=40),
        showlegend=True,
        legend=dict(orientation="v", yanchor="middle", y=0.5, xanchor="left", x=1.02)
    )
    
    fig.write_image("borghi_empty.pdf")
    fig.show()

# Load and filter data
df = pd.read_csv('/Users/raja/Library/CloudStorage/OneDrive-SharedLibraries-McGillUniversity/CombustionDynamicsF2025 - General/04. Code - git repo/src/test_data.csv')
df = df[df['dP_percent'] < 0]
df = df[df['M'] < 0.006]
df = df[df['BR'] < 0.5]
df = df[df['U0'] == 90]
df = df[df['PHI'] == 0.47]
df = df[df['L_premixer'] == 0.10674]


# Create combined plot
#create_combined_plots(df, xScatter="lt_Lf", yScatter="ui_SL")


import plotly.graph_objects as go
import plotly.express as px

def BorghiPlot(df: pd.DataFrame, xScatter: str, yScatter: str) -> None:

    fig = go.Figure()

    # GRAPHS
    # laminar curve
    x_lam = np.array([0.1, 10])
    y_lam = np.array([10, 0.1])
    fig.add_trace(go.Scatter(x=x_lam, y=y_lam, mode="lines",
        line=dict(color="black", width=1),
        name=r"$Re_t = 1$"))
    fig.add_annotation(x=np.log10(0.5), y=np.log10(0.3), text="Laminar",
        showarrow=False, font=dict(size=14))
    
    # u_prime/S_L = 1
    x_2 = np.array([1.0, 1000])
    y_2 = np.array([1.0, 1])
    fig.add_trace(go.Scatter(x=x_2, y=y_2, mode="lines",
        line=dict(color="black", width=1),
        name=r"$u'/S_L = 1$"))
    fig.add_annotation(x=np.log10(400), y=np.log10(0.85), text=r"$u'/S_L = 1$",
        showarrow=False, font=dict(size=14))


    # Now we make plot for Ka = 1
    x_3 = np.linspace(1.0, 1000, 1000)
    y_3 = x_3**(1/3)
    fig.add_trace(go.Scatter(x=x_3, y=y_3, mode="lines",
        line=dict(color="black", width=1),
        name="Ka = 1"))
    fig.add_annotation(x=np.log10(400), y=np.log10(10), text="Ka = 1",
        showarrow=False, font=dict(size=14))
    
    # Now we make plot for Ka = 10
    x_3a = np.linspace(0.317, 1000, 1000)
    y_3a = (10 * (x_3a)**(1/2))**(2/3)
    fig.add_trace(go.Scatter(x=x_3a, y=y_3a, mode="lines",
        line=dict(color="black", width=1, dash="dash"),
        name="Ka = 10"))
    fig.add_annotation(x=np.log10(2), y=np.log10(36), text="Ka = 100",
        showarrow=False, font=dict(size=14))


    # Now we make plot for Ka = 100
    x_3b = np.linspace(0.1, 1000, 1000)
    y_3b = (100 * (x_3b)**(1/2))**(2/3)
    fig.add_trace(go.Scatter(x=x_3b, y=y_3b, mode="lines",
        line=dict(color="black", width=1),
        name="Ka = 100"))
    fig.add_annotation(x=np.log10(2), y=np.log10(36), text="Ka = 100",
        showarrow=False, font=dict(size=14))


    # Now we make plot for Da = 1
    x_4 = np.linspace(1.0, 1000, 1000)
    y_4 = x_4
    fig.add_trace(go.Scatter(x=x_4, y=y_4, mode="lines",
        line=dict(color="black", width=1),
        name="Da = 1"))
    fig.add_annotation(x=np.log10(90), y=np.log10(150), text="Da = 1",
        showarrow=False, font=dict(size=14))
    
    df['hover_text'] = df.apply(lambda row: f"U0: {row['U0']}<br>PHI: {row['PHI']}", axis=1)

    # Scatter data
    scatter_fig = px.scatter(
        df, x=xScatter, y=yScatter, color="Stabilization", hover_data=["BR"])
        #size='U0')
    for trace in scatter_fig.data:
        fig.add_trace(trace)

    # set the symbol as star
    #fig.update_traces(marker=dict(color="red", size=16, symbol='star'))

    #hide_types = ["Bluff Body Flame", "Swirl Flame"]     # change as needed
    #for trace in fig.data:
    #    if trace.name in hide_types:
    #        trace.visible = "legendonly"

    #fig.update_traces(marker=dict(size=12))
    fig.update_xaxes(
        type="log",
        autorange=False, range=[-1, 3], 
        showgrid=False,showline=True,mirror=True,automargin=True,zeroline=False,
        title=r"$l_t/\delta_L$"
        )
    fig.update_yaxes(
        type="log",
        autorange=False, range=[-1, 2.8],
        showgrid=False,showline=True, mirror=True,automargin=True,zeroline=False ,
        title=r"$u'/S_L$"
        )

    fig.update_layout(
        title=dict(text="Borghi Diagram", x=0.5, y=0.98, font=dict(size=16)),
        #legend=dict(orientation="h", yanchor="bottom", xanchor="right", x=1),
        margin=dict(l=20, r=20, t=40, b=40),
        template="simple_white",
        width=800, height=600,
        showlegend=True
    )
    fig.write_image("BP_empty.pdf")
    fig.show()

#df = pd.read_csv('/Users/raja/Library/CloudStorage/OneDrive-SharedLibraries-McGillUniversity/CombustionDynamicsF2025 - General/04. Code - git repo/src/test_data.csv')
#df = df[df['dP_percent'] < 1.6]
#df = df[df['M'] < 0.006]
#df = df[df['BR'] < 0.5]
#df = df[df['U0'] == 90]
#df = df[df['PHI'] == 0.47]
BorghiPlot(df, xScatter="lt_Lf", yScatter="ui_SL")
