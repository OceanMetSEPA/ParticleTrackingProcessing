import re

def parse_xml_row(line):
    """ Parses an XML row and extracts the variable name and value. """
    line = line.strip()
    variable_name = ''
    val = None
    
    try:
        if '=' in line:
            # Handle key-value pairs like "Key=Value"
            parts = line.split('=')
            variable_name = parts[0].strip('<>')
            val = parts[1].strip('<>/')
        else:
            # Handle XML tag-value pairs like <DateTime>2022-11-13 00:30:00</DateTime>
            match = re.search(r'<(.*?)>(.*?)</\1>', line)
            if match:
                variable_name, val = match.groups()
        
        # Convert numerical values
        try:
            val=val.replace('"','')
            val = float(val)
        except ValueError:
            pass
        
    except Exception:# as e:
#        print(f"Error parsing line: {line} -> {e}")
         'oh well'
    
    return variable_name, val