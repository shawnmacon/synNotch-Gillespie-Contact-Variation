"""
Simulation Database Configuration and Accessing
Shawn Macon
09.09.2025
"""

import sqlite3


def create_schema():
    conn = sqlite3.connect('simulation_data.db')
    cursor = conn.cursor()

    # Create table for simulation parameters (now includes redraw_interval)
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS parameters (
        synthesis_constant REAL NOT NULL,
        decay_constant REAL NOT NULL,
        LoP_average REAL NOT NULL,
        LoP_deviation_value REAL NOT NULL,
        redraw_interval REAL NOT NULL,
        PRIMARY KEY (synthesis_constant, decay_constant, LoP_average, LoP_deviation_value, redraw_interval)
    );
    ''')

    # Create table for GFP concentration data (with redraw_interval and trial_id)
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS gfp_concentration_data (
        synthesis_constant REAL,
        decay_constant REAL,
        LoP_average REAL,
        LoP_deviation_value REAL,
        redraw_interval REAL,
        trial_id INTEGER NOT NULL,
        time_point REAL,
        gfp_concentration REAL,
        PRIMARY KEY (synthesis_constant, decay_constant, LoP_average, LoP_deviation_value, redraw_interval, trial_id, time_point),
        FOREIGN KEY (synthesis_constant, decay_constant, LoP_average, LoP_deviation_value, redraw_interval)
            REFERENCES parameters (synthesis_constant, decay_constant, LoP_average, LoP_deviation_value, redraw_interval)
    );
    ''')

    conn.commit()
    conn.close()


def insert_simulation_data(parameters, gfp_data, trial_id):
    """
    parameters = (synthesis_constant, decay_constant, LoP_average, LoP_deviation_value, redraw_interval)
    gfp_data = list of (time_point, gfp_concentration)
    """
    conn = sqlite3.connect('simulation_data.db')
    cursor = conn.cursor()
    cursor.execute("BEGIN TRANSACTION;")

    try:
        # Insert simulation parameters
        cursor.execute(''' 
        INSERT OR REPLACE INTO parameters (synthesis_constant, decay_constant, LoP_average, LoP_deviation_value, redraw_interval)
        VALUES (?, ?, ?, ?, ?)
        ''', parameters)

        # Prepare GFP data for bulk insert
        gfp_data_to_insert = [
            (parameters[0], parameters[1], parameters[2], parameters[3], parameters[4],
             trial_id, time_point, gfp_concentration)
            for time_point, gfp_concentration in gfp_data
        ]

        # Insert GFP concentration data
        cursor.executemany(''' 
        INSERT OR REPLACE INTO gfp_concentration_data 
        (synthesis_constant, decay_constant, LoP_average, LoP_deviation_value, redraw_interval, trial_id, time_point, gfp_concentration)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        ''', gfp_data_to_insert)

        conn.commit()

    except Exception as e:
        print(f"Error inserting data: {e}")
        conn.rollback()
    finally:
        conn.close()


def insert_multiple_simulations(simulations):
    """
    Inserts multiple simulations into the database in one batch.
    simulations = list of (parameters, gfp_data)
    where parameters = (synthesis_constant, decay_constant, LoP_average, LoP_deviation_value, redraw_interval)
    and gfp_data = list of (trial_id, [(time_point, gfp_concentration), ...])
    """
    conn = sqlite3.connect('simulation_data.db')
    cursor = conn.cursor()
    cursor.execute('BEGIN TRANSACTION;')

    try:
        for parameters, trials in simulations:
            # Insert parameters once
            cursor.execute('''
            INSERT OR REPLACE INTO parameters 
            (synthesis_constant, decay_constant, LoP_average, LoP_deviation_value, redraw_interval)
            VALUES (?, ?, ?, ?, ?)
            ''', parameters)

            # Each element in `trials` is (trial_id, [(t, gfp)])
            for trial_id, gfp_data in trials:
                gfp_data_to_insert = [
                    (parameters[0], parameters[1], parameters[2], parameters[3], parameters[4],
                     trial_id, time_point, gfp_concentration)
                    for time_point, gfp_concentration in gfp_data
                ]

                cursor.executemany('''
                INSERT OR REPLACE INTO gfp_concentration_data 
                (synthesis_constant, decay_constant, LoP_average, LoP_deviation_value, redraw_interval, trial_id, time_point, gfp_concentration)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                ''', gfp_data_to_insert)

        conn.commit()

    except Exception as e:
        print(f"Error inserting multiple simulations: {e}")
        conn.rollback()
    finally:
        conn.close()


def get_simulation_results(synthesis_constant, decay_constant, LoP_average, LoP_deviation_value, redraw_interval, trial_id=None):
    conn = sqlite3.connect('simulation_data.db')
    cursor = conn.cursor()

    if trial_id is not None:
        cursor.execute('''
        SELECT time_point, gfp_concentration
        FROM gfp_concentration_data
        WHERE synthesis_constant = ?
        AND decay_constant = ?
        AND LoP_average = ?
        AND LoP_deviation_value = ?
        AND redraw_interval = ?
        AND trial_id = ?
        ORDER BY time_point
        ''', (synthesis_constant, decay_constant, LoP_average, LoP_deviation_value, redraw_interval, trial_id))
        results = cursor.fetchall()
        conn.close()
        return [(t, g) for t, g in results]

    else:
        cursor.execute('''
        SELECT trial_id, time_point, gfp_concentration
        FROM gfp_concentration_data
        WHERE synthesis_constant = ?
        AND decay_constant = ?
        AND LoP_average = ?
        AND LoP_deviation_value = ?
        AND redraw_interval = ?
        ORDER BY trial_id, time_point
        ''', (synthesis_constant, decay_constant, LoP_average, LoP_deviation_value, redraw_interval))
        results = cursor.fetchall()

        trials_data = {}
        for trial_id, time_point, gfp_concentration in results:
            trials_data.setdefault(trial_id, []).append((time_point, gfp_concentration))

        conn.close()
        return trials_data