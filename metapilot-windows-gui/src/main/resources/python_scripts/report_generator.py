import os
import json
from typing import Dict, List, Optional

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from flask import Flask, jsonify, render_template, request
from flask_cors import CORS

app = Flask(__name__)
CORS(app)

class ReportGenerator:
    def __init__(self, output_dir: str = "reports"):
        """Initialize the report generator.
        
        Args:
            output_dir: Directory to save generated reports
        """
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)

    def create_correlation_plot(self, data: pd.DataFrame, output_file: str) -> str:
        """Create a correlation heatmap plot.
        
        Args:
            data: DataFrame containing numerical data
            output_file: Name of the output file
            
        Returns:
            Path to the generated plot
        """
        corr_matrix = data.corr()
        fig = px.imshow(
            corr_matrix,
            title="Correlation Matrix",
            color_continuous_scale="RdBu_r"
        )
        
        output_path = os.path.join(self.output_dir, output_file)
        fig.write_html(output_path)
        return output_path

    def create_distribution_plot(
        self,
        data: pd.DataFrame,
        column: str,
        output_file: str,
        bins: int = 30
    ) -> str:
        """Create a distribution plot for a specific column.
        
        Args:
            data: DataFrame containing the data
            column: Column name to plot
            output_file: Name of the output file
            bins: Number of histogram bins
            
        Returns:
            Path to the generated plot
        """
        fig = px.histogram(
            data,
            x=column,
            nbins=bins,
            title=f"Distribution of {column}"
        )
        
        output_path = os.path.join(self.output_dir, output_file)
        fig.write_html(output_path)
        return output_path

    def generate_report(
        self,
        data: pd.DataFrame,
        report_name: str,
        plots: Optional[List[Dict]] = None
    ) -> Dict:
        """Generate a complete report with multiple visualizations.
        
        Args:
            data: DataFrame containing the data
            report_name: Name of the report
            plots: List of plot configurations
            
        Returns:
            Dictionary containing report metadata and plot paths
        """
        if plots is None:
            plots = [
                {"type": "correlation", "output": f"{report_name}_correlation.html"},
                {"type": "distribution", "column": data.columns[0], "output": f"{report_name}_distribution.html"}
            ]
        
        report_data = {
            "name": report_name,
            "plots": []
        }
        
        for plot_config in plots:
            if plot_config["type"] == "correlation":
                plot_path = self.create_correlation_plot(
                    data,
                    plot_config["output"]
                )
            elif plot_config["type"] == "distribution":
                plot_path = self.create_distribution_plot(
                    data,
                    plot_config["column"],
                    plot_config["output"]
                )
            
            report_data["plots"].append({
                "type": plot_config["type"],
                "path": plot_path
            })
        
        return report_data

@app.route("/generate_report", methods=["POST"])
def generate_report_endpoint():
    """API endpoint for generating reports."""
    try:
        data = request.json
        df = pd.DataFrame(data["data"])
        report_name = data.get("report_name", "default_report")
        plots = data.get("plots")
        
        generator = ReportGenerator()
        report = generator.generate_report(df, report_name, plots)
        
        return jsonify({"status": "success", "report": report})
    except Exception as e:
        return jsonify({"status": "error", "message": str(e)}), 500

if __name__ == "__main__":
    app.run(debug=True, port=5000) 