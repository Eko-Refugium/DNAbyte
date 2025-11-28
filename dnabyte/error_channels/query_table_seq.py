import psycopg2
import pandas as pd

def run_query():
    # Database connection parameters
    db_params = {
        'dbname': 'dna_sim',
        'user': 'dna_sim',
        'password': 'test1337'
    }

    # Connect to the PostgreSQL database
    conn = psycopg2.connect(**db_params)
    cursor = conn.cursor()

    # Define the query to fetch the entire table
    query = "SELECT * FROM public.seq_err_rates"

    try:
        # Execute the query
        cursor.execute(query)
        result = cursor.fetchall()

        # Fetch column names
        colnames = [desc[0] for desc in cursor.description]

        # Create a DataFrame
        df = pd.DataFrame(result, columns=colnames)

        # Save the DataFrame as a pickle file
        df.to_pickle("seq_table.pkl")

        # Save the DataFrame as a JSON file
        df.to_json("seq_table.json", orient='records', lines=True)


    except Exception as e:
        #logger.info("An error occurred:", e)
        pass
    
    finally:
        # Close the cursor and connection
        cursor.close()
        conn.close()

if __name__ == "__main__":
    run_query()