import sqlite3
import pandas as pd

DB_PATH = "genome.db"
CSV_PATH = "data/raw/sample_genome_metrics.csv"

def load_data():
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    df = pd.read_csv(CSV_PATH)

    # Convert wide → long (key-value format)
    melted_df = df.melt(
        id_vars=["assembly_id"],
        var_name="metric_name",
        value_name="value"
    )

    # Insert into DB
    for _, row in melted_df.iterrows():
        cursor.execute("""
            INSERT INTO assembly_metrics (assembly_id, metric_name, value)
            VALUES (?, ?, ?)
        """, (row["assembly_id"], row["metric_name"], row["value"]))

    conn.commit()
    conn.close()

if __name__ == "__main__":
    load_data()