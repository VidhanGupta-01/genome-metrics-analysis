import sqlite3
import pandas as pd

DB_PATH = "genome.db"

class TransformationService:

    def fetch_metrics(self):
        conn = sqlite3.connect(DB_PATH)

        query = """
        SELECT assembly_id, metric_name, value
        FROM assembly_metrics
        """

        df = pd.read_sql(query, conn)
        conn.close()
        return df

    def pivot_metrics(self, df):
        df_wide = df.pivot(
            index="assembly_id",
            columns="metric_name",
            values="value"
        )
        return df_wide

    def clean_data(self, df):
        df = df.fillna(df.median(numeric_only=True))
        return df

    def transform(self):
        df = self.fetch_metrics()
        df = self.pivot_metrics(df)
        df = self.clean_data(df)
        return df