import numpy as np
import pandas as pd
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
        name="Re_t = 1"))
    fig.add_annotation(x=np.log10(0.5), y=np.log10(0.3), text="Laminar",
        showarrow=False, font=dict(size=14))
    
    # u_prime/S_L = 1
    x_2 = np.array([1.0, 1000])
    y_2 = np.array([1.0, 1])
    fig.add_trace(go.Scatter(x=x_2, y=y_2, mode="lines",
        line=dict(color="black", width=1),
        name="Re_t = 1"))
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
    fig.add_annotation(x=np.log10(90), y=np.log10(50), text="Da = 1",
        showarrow=False, font=dict(size=14))
    

    # Scatter data
    scatter_fig = px.scatter(
        df, x=xScatter, y=yScatter, color="U0",)
        #size='U0')
    for trace in scatter_fig.data:
        fig.add_trace(trace)


    fig.update_xaxes(
        type="log",
        autorange=False, range=[-1, 3], 
        showgrid=False,showline=True,mirror=True,automargin=True,zeroline=False,
        title=r"$l_t/\delta_L$"
        )
    fig.update_yaxes(
        type="log",
        autorange=False, range=[-1, 2],
        showgrid=False,showline=True, mirror=True,automargin=True,zeroline=False ,
        title=r"$u'/S_L$"
        )

    fig.update_layout(
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        margin=dict(l=60, r=20, t=80, b=60),
        template="simple_white",
        width=840, height=640,
        showlegend=False
    )
    fig.write_html("BP.html")
    fig.show()

def test_borghi():
    df = pd.read_csv('test_data.csv')
    BorghiPlot(df, 'U0', 'S_L')

if __name__ == "__main__":
    test_borghi()