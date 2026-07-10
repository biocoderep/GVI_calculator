import re
newick = "(OP476775.1|06-Sep-2009:0.0112056488)"
print(re.sub(r'([A-Za-z0-9_\.]+)\|[^:,()]+', r'\1', newick))
