import sqlite3
import json
import logging
from config import Config
import os

logger = logging.getLogger(__name__)

DB_PATH = os.path.join(Config.BASE_DIR, 'jobs.db')

def get_db():
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn

def init_db():
    try:
        with get_db() as conn:
            conn.execute('''
                CREATE TABLE IF NOT EXISTS jobs (
                    id TEXT PRIMARY KEY,
                    status TEXT NOT NULL,
                    progress INTEGER NOT NULL DEFAULT 0,
                    current_file TEXT,
                    current_step TEXT,
                    results TEXT,
                    message TEXT,
                    csv_path TEXT
                )
            ''')
            conn.commit()
    except Exception as e:
        logger.error(f"Error initializing jobs database: {e}")

def create_job(job_id):
    try:
        with get_db() as conn:
            conn.execute(
                'INSERT INTO jobs (id, status, progress, current_file, current_step) VALUES (?, ?, ?, ?, ?)',
                (job_id, 'queued', 0, 'Starting...', 'Initializing Engine...')
            )
            conn.commit()
    except Exception as e:
        logger.error(f"Error creating job {job_id}: {e}")

def update_job(job_id, status=None, progress=None, current_file=None, current_step=None, results=None, message=None, csv_path=None):
    updates = []
    params = []
    
    if status is not None:
        updates.append("status = ?")
        params.append(status)
    if progress is not None:
        updates.append("progress = ?")
        params.append(progress)
    if current_file is not None:
        updates.append("current_file = ?")
        params.append(current_file)
    if current_step is not None:
        updates.append("current_step = ?")
        params.append(current_step)
    if results is not None:
        updates.append("results = ?")
        params.append(json.dumps(results))
    if message is not None:
        updates.append("message = ?")
        params.append(message)
    if csv_path is not None:
        updates.append("csv_path = ?")
        params.append(csv_path)
        
    if not updates:
        return
        
    params.append(job_id)
    query = f"UPDATE jobs SET {', '.join(updates)} WHERE id = ?"
    
    try:
        with get_db() as conn:
            conn.execute(query, params)
            conn.commit()
    except Exception as e:
        logger.error(f"Error updating job {job_id}: {e}")

def get_job(job_id):
    try:
        with get_db() as conn:
            row = conn.execute('SELECT * FROM jobs WHERE id = ?', (job_id,)).fetchone()
            if row:
                job_dict = dict(row)
                if job_dict['results']:
                    job_dict['results'] = json.loads(job_dict['results'])
                else:
                    job_dict['results'] = []
                return job_dict
            return None
    except Exception as e:
        logger.error(f"Error retrieving job {job_id}: {e}")
        return None
