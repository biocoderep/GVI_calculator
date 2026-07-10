import os

with open("E:/Inlfuenza/avian_influenza/GVI_index_calculator/scientific_engine/external_tools.py", "r") as f:
    content = f.read()

ess_func = """
def calculate_ess(x):
    import numpy as np
    n = len(x)
    if n < 2: return 0.0
    mu = np.mean(x)
    var = np.var(x)
    if var == 0: return float(n)
    auto_cov = []
    for lag in range(1, min(n//2, 1000)):
        cov = np.sum((x[:-lag] - mu) * (x[lag:] - mu)) / n
        r = cov / var
        if r < 0.05: break
        auto_cov.append(r)
    tau = 1 + 2 * sum(auto_cov)
    return max(1.0, n / tau)

"""

content = content.replace("logger = logging.getLogger(__name__)\n", "logger = logging.getLogger(__name__)\n" + ess_func)

beast_old = """        cmd = [beast_bin, "-overwrite", xml_path]
        if Config.USE_WSL_FALLBACK: cmd = ["wsl"] + cmd
            
        try:
            logger.info(f"Starting BEAST MCMC on chunk {i+1}/{len(chunk_paths)}")
            work_dir = os.path.dirname(chunk)
            subprocess.run(cmd, cwd=work_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=None, text=True)
            
            log_path = os.path.join(work_dir, log_file)
            if os.path.exists(log_path):
                import pandas as pd
                df = pd.read_csv(log_path, sep='\\t', comment='#')
                burn_in = int(len(df) * 0.1)
                if 'recombinationRate' in df.columns:
                    valid_samples = df['recombinationRate'].iloc[burn_in:]
                    if not valid_samples.empty: chunk_recomb = float(valid_samples.mean())
                if 'ucldMean.c:alignment' in df.columns:
                    valid_samples = df['ucldMean.c:alignment'].iloc[burn_in:]
                    if not valid_samples.empty: chunk_evo = float(valid_samples.mean())
                    
                recombination_rates.append(chunk_recomb)
                beast_evo_rates.append(chunk_evo)
                weights.append(seq_count)
        except Exception as e:
            logger.error(f"BEAST chunk {i+1} failed: {e}")"""

beast_new = """        cmd = [beast_bin, "-overwrite", xml_path]
        if Config.USE_WSL_FALLBACK: cmd = ["wsl"] + cmd
            
        loop_count = 0
        max_loops = 5
        ess_achieved = False
        
        while loop_count < max_loops and not ess_achieved:
            try:
                logger.info(f"Starting BEAST MCMC on chunk {i+1}/{len(chunk_paths)} (Loop {loop_count+1})")
                work_dir = os.path.dirname(chunk)
                subprocess.run(cmd, cwd=work_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=None, text=True)
                
                log_path = os.path.join(work_dir, log_file)
                if os.path.exists(log_path):
                    import pandas as pd
                    df = pd.read_csv(log_path, sep='\\t', comment='#')
                    burn_in = int(len(df) * 0.1)
                    
                    recomb_ess = 200.0
                    evo_ess = 200.0
                    
                    if 'recombinationRate' in df.columns:
                        valid_samples = df['recombinationRate'].iloc[burn_in:]
                        if not valid_samples.empty:
                            chunk_recomb = float(valid_samples.mean())
                            recomb_ess = calculate_ess(valid_samples.values)
                            
                    if 'ucldMean.c:alignment' in df.columns:
                        valid_samples = df['ucldMean.c:alignment'].iloc[burn_in:]
                        if not valid_samples.empty:
                            chunk_evo = float(valid_samples.mean())
                            evo_ess = calculate_ess(valid_samples.values)
                            
                    logger.info(f"Chunk {i+1} ESS - Recomb: {recomb_ess:.1f}, Evo: {evo_ess:.1f}")
                    if recomb_ess < 200 or evo_ess < 200:
                        logger.warning(f"ESS < 200. Resuming BEAST MCMC...")
                        cmd = [beast_bin, "-resume", xml_path]
                        if Config.USE_WSL_FALLBACK: cmd = ["wsl"] + cmd
                        loop_count += 1
                    else:
                        ess_achieved = True
                        recombination_rates.append(chunk_recomb)
                        beast_evo_rates.append(chunk_evo)
                        weights.append(seq_count)
                        break
            except Exception as e:
                logger.error(f"BEAST chunk {i+1} failed: {e}")
                break
                
        if not ess_achieved and loop_count >= max_loops:
            logger.warning(f"Max loops reached for chunk {i+1}. Using sub-optimal ESS.")
            recombination_rates.append(chunk_recomb)
            beast_evo_rates.append(chunk_evo)
            weights.append(seq_count)"""

content = content.replace(beast_old, beast_new)

bdsky_old = """        cmd = [beast_bin, "-overwrite", xml_path]
        if Config.USE_WSL_FALLBACK: cmd = ["wsl"] + cmd
            
        try:
            logger.info(f"Starting BDSKY MCMC on chunk {i+1}/{len(chunk_paths)}")
            work_dir = os.path.dirname(chunk)
            subprocess.run(cmd, cwd=work_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=None, text=True)
            
            log_path = os.path.join(work_dir, log_file)
            if os.path.exists(log_path):
                import pandas as pd
                df = pd.read_csv(log_path, sep='\\t', comment='#')
                if 'reproductiveNumber' in df.columns:
                    burn_in = int(len(df) * 0.1)
                    valid_samples = df['reproductiveNumber'].iloc[burn_in:]
                    if not valid_samples.empty:
                        chunk_re = float(valid_samples.mean())
                        re_rates.append(chunk_re)
                        weights.append(seq_count)
        except Exception as e:
            pass"""

bdsky_new = """        cmd = [beast_bin, "-overwrite", xml_path]
        if Config.USE_WSL_FALLBACK: cmd = ["wsl"] + cmd
            
        loop_count = 0
        max_loops = 5
        ess_achieved = False
        
        while loop_count < max_loops and not ess_achieved:
            try:
                logger.info(f"Starting BDSKY MCMC on chunk {i+1}/{len(chunk_paths)} (Loop {loop_count+1})")
                work_dir = os.path.dirname(chunk)
                subprocess.run(cmd, cwd=work_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=None, text=True)
                
                log_path = os.path.join(work_dir, log_file)
                if os.path.exists(log_path):
                    import pandas as pd
                    df = pd.read_csv(log_path, sep='\\t', comment='#')
                    re_ess = 200.0
                    if 'reproductiveNumber' in df.columns:
                        burn_in = int(len(df) * 0.1)
                        valid_samples = df['reproductiveNumber'].iloc[burn_in:]
                        if not valid_samples.empty:
                            chunk_re = float(valid_samples.mean())
                            re_ess = calculate_ess(valid_samples.values)
                            
                    logger.info(f"Chunk {i+1} ESS - Re: {re_ess:.1f}")
                    if re_ess < 200:
                        logger.warning(f"ESS < 200. Resuming BDSKY MCMC...")
                        cmd = [beast_bin, "-resume", xml_path]
                        if Config.USE_WSL_FALLBACK: cmd = ["wsl"] + cmd
                        loop_count += 1
                    else:
                        ess_achieved = True
                        re_rates.append(chunk_re)
                        weights.append(seq_count)
                        break
            except Exception as e:
                break
                
        if not ess_achieved and loop_count >= max_loops:
            logger.warning(f"Max loops reached for BDSKY chunk {i+1}. Using sub-optimal ESS.")
            re_rates.append(chunk_re)
            weights.append(seq_count)"""

content = content.replace(bdsky_old, bdsky_new)

with open("E:/Inlfuenza/avian_influenza/GVI_index_calculator/scientific_engine/external_tools.py", "w") as f:
    f.write(content)
