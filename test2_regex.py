import re
newick = "(PQ768018|India__Goa|cattle|2013:0.111, PQ768017|India__Jharkhand|cattle|2013:0.222);"
print(re.sub(r'([A-Za-z0-9_\.]+)\|[^:,()]+', r'\1', newick))
