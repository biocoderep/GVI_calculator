import os

class Config:
    # Flask settings
    SECRET_KEY = os.environ.get('SECRET_KEY', 'dev-secret-key')
    
    # Paths
    BASE_DIR = os.path.abspath(os.path.dirname(__file__))
    UPLOAD_FOLDER = os.path.join(BASE_DIR, 'uploads')
    TREES_FOLDER = os.path.join(BASE_DIR, 'static', 'trees')
    
    # Upload limits
    MAX_CONTENT_LENGTH = 50 * 1024 * 1024  # 50MB max upload
    
    # External Tools
    IQTREE_BIN = "iqtree2"
    BEAST_BIN = "beast"
    PHIPACK_BIN = "Phi"
    
    import platform
    # If running in Docker, we are already natively on Linux, so no WSL is needed.
    if os.environ.get('DOCKER_ENV') == 'true':
        USE_WSL_FALLBACK = False
    else:
        USE_WSL_FALLBACK = (platform.system() == "Windows")
