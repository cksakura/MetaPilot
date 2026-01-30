"""
Create word cloud visualization from word_cloud_data.csv
"""
import pandas as pd
import matplotlib.pyplot as plt
from wordcloud import WordCloud
import sys
import os

def create_wordcloud(csv_file, output_file='wordcloud.png'):
    """Create word cloud from CSV file"""
    try:
        # Read word count data
        df = pd.read_csv(csv_file)
        
        # Create dictionary from dataframe
        word_dict = dict(zip(df['word'], df['count']))
        
        # Create word cloud
        wordcloud = WordCloud(
            width=1200,
            height=800,
            background_color='white',
            max_words=200,
            colormap='viridis',
            relative_scaling=0.5,
            min_font_size=10
        ).generate_from_frequencies(word_dict)
        
        # Plot
        plt.figure(figsize=(15, 10))
        plt.imshow(wordcloud, interpolation='bilinear')
        plt.axis('off')
        plt.title('PRIDE Dataset Word Cloud', fontsize=20, pad=20)
        plt.tight_layout(pad=0)
        
        # Save
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Word cloud saved to: {output_file}")
        
        # Also create HTML version for interactive viewing
        html_file = output_file.replace('.png', '.html')
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>PRIDE Dataset Word Cloud</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            text-align: center;
            background-color: #f5f5f5;
        }}
        h1 {{
            color: #333;
            margin: 20px;
        }}
        img {{
            max-width: 100%;
            height: auto;
            box-shadow: 0 4px 8px rgba(0,0,0,0.1);
        }}
    </style>
</head>
<body>
    <h1>PRIDE Dataset Word Cloud</h1>
    <img src="{os.path.basename(output_file)}" alt="Word Cloud">
</body>
</html>
"""
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        print(f"HTML file saved to: {html_file}")
        
    except Exception as e:
        print(f"Error creating word cloud: {e}")
        print("Make sure you have installed: pip install wordcloud matplotlib pandas")
        sys.exit(1)

if __name__ == "__main__":
    csv_file = "word_cloud_data.csv"
    if len(sys.argv) > 1:
        csv_file = sys.argv[1]
    
    if not os.path.exists(csv_file):
        print(f"Error: {csv_file} not found")
        print("Run AnalyzePrideDatasets.java first to generate the CSV file")
        sys.exit(1)
    
    output_file = "wordcloud.png"
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    
    create_wordcloud(csv_file, output_file)

