import sys
content = open('scientific_engine/external_tools.py').read()
content = content.replace('beast_bin = "/mnt/e/Inlfuenza/avian_influenza/GVI_index_calculator/tools/beast/bin/beast"', 'beast_bin = Config.BEAST_BIN')
open('scientific_engine/external_tools.py', 'w').write(content)
