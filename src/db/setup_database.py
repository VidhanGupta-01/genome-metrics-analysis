import sqlite3

def create_database():
    conn = sqlite3.connect("genome.db")
    cursor = conn.cursor()

    cursor.execute("""
    CREATE TABLE IF NOT EXISTS assembly_metrics (
        assembly_id TEXT,
        metric_name TEXT,
        value REAL
    )
    """)

    conn.commit()
    conn.close()

if __name__ == "__main__":
    create_database()