import psycopg2
import pandas as pd

def run_query():
    # Database connection parameters
    db_params = {
        'dbname': 'dna_sim',
        'user': 'fabian'
    }

    # Connect to the PostgreSQL database
    conn = psycopg2.connect(**db_params)
    cursor = conn.cursor()

    # Define the query to fetch the entire table
    query = "SELECT * FROM synth_err_rates"

    try:
        # Execute the query
        cursor.execute(query)
        result = cursor.fetchall()

        # Fetch column names
        colnames = [desc[0] for desc in cursor.description]

        # Create a DataFrame
        df = pd.DataFrame(result, columns=colnames)

        # Save the DataFrame as a pickle file
        df.to_pickle("syn_table.pkl")

    except Exception as e:
        # logger.info(f"Error: {e}")
        pass

    finally:
        # Close the cursor and connection
        cursor.close()
        conn.close()

if __name__ == "__main__":
    run_query()