from flask import Flask, render_template, request, jsonify
import os
import uuid
import threading
import logging
from werkzeug.utils import secure_filename
from bio_engine import process_fasta_file
import pandas as pd
from config import Config
import job_manager

# Setup structured logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

app = Flask(__name__)
app.config.from_object(Config)

os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
os.makedirs(app.config['TREES_FOLDER'], exist_ok=True)

# Initialize Jobs DB
job_manager.init_db()

@app.route('/')
def index():
    return render_template('index.html')

def background_processor(job_id, file_paths, custom_weights=None):
    logger.info(f"Starting background job {job_id} with {len(file_paths)} files")
    job_manager.update_job(job_id, status='running')
    results = []
    
    total = len(file_paths)
    for i, filepath in enumerate(file_paths):
        current_file = os.path.basename(filepath)
        job_manager.update_job(job_id, current_file=current_file)
        try:
            res = process_fasta_file(filepath, custom_weights, job_id)
            if res:
                results.append(res)
        except Exception as e:
            logger.error(f"Error processing {filepath}: {str(e)}", exc_info=True)
        
        progress = int(((i + 1) / total) * 100)
        job_manager.update_job(job_id, progress=progress)
    
    # Save results to CSV
    if results:
        df = pd.DataFrame(results)
        # Sort by year if possible
        try:
            df['Year_Int'] = pd.to_numeric(df['Year'], errors='coerce')
            df = df.sort_values(by='Year_Int').drop(columns=['Year_Int'])
        except:
            pass
            
        csv_path = os.path.join(app.config['UPLOAD_FOLDER'], f"{job_id}_results.csv")
        df.to_csv(csv_path, index=False)
        job_manager.update_job(
            job_id, 
            status='completed', 
            results=df.to_dict(orient='records'),
            csv_path=csv_path
        )
        logger.info(f"Job {job_id} completed successfully. Saved to {csv_path}")
    else:
        job_manager.update_job(
            job_id,
            status='error',
            message='No valid sequences processed.'
        )
        logger.warning(f"Job {job_id} failed: No valid sequences processed.")

@app.route('/upload', methods=['POST'])
def upload_files():
    logger.info("Received upload request")
    if 'files[]' not in request.files:
        logger.warning("Upload failed: No files[] in request")
        return jsonify({'error': 'No files uploaded'}), 400
        
    files = request.files.getlist('files[]')
    if not files or files[0].filename == '':
        logger.warning("Upload failed: No selected files")
        return jsonify({'error': 'No selected files'}), 400
        
    job_id = str(uuid.uuid4())
    job_manager.create_job(job_id)
    
    saved_paths = []
    for file in files:
        if file and file.filename.endswith(('.fasta', '.fa')):
            filename = secure_filename(file.filename)
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], f"{job_id}_{filename}")
            file.save(filepath)
            saved_paths.append(filepath)
            
    if not saved_paths:
        logger.warning(f"Upload failed: No valid FASTA files found for job {job_id}")
        return jsonify({'error': 'No valid FASTA files found'}), 400
        
    logger.info(f"Created job {job_id} with {len(saved_paths)} files")
    
    custom_weights_str = request.form.get('weights')
    custom_weights = {}
    if custom_weights_str:
        import json
        try:
            custom_weights = json.loads(custom_weights_str)
        except Exception:
            pass

    # Start background processing thread
    thread = threading.Thread(target=background_processor, args=(job_id, saved_paths, custom_weights))
    thread.start()
    
    return jsonify({'job_id': job_id})

@app.route('/status/<job_id>', methods=['GET'])
def check_status(job_id):
    job = job_manager.get_job(job_id)
    if not job:
        return jsonify({'error': 'Job not found'}), 404
    return jsonify(job)

if __name__ == '__main__':
    app.run(host='0.0.0.0', debug=True, port=5000)
