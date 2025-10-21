#!/usr/bin/env python3
import dash
from dash import html, dcc
import dash_bootstrap_components as dbc

# Create simple test app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

app.layout = html.Div([
    html.H1("GeoCLUSTER Test App"),
    html.P("If you can see this, Dash is working correctly!"),
    dbc.Button("Test Button", color="primary"),
])

if __name__ == "__main__":
    print("Starting simple test app on http://localhost:3000")
    app.run_server(debug=True, port=3000, host='127.0.0.1')
