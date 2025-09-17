# run_xml_processing.py
import sys
from pathlib import Path

# Ensure xmltools is on Python's path
xmltools_path = r"C:\Python\xmlStuff"  # folder containing xmltools
if xmltools_path not in sys.path:
    sys.path.insert(0, xmltools_path)

from xmltools.process_xml_folder import process_xml_folder

#%%
%%time

# Folder containing XML files (update as needed)
model_path = r"\\cloudnet.sepa.org.uk\sepadata\FRMData\Marine\SeaLiceScreening"
model_path = Path(model_path)

if not model_path.is_dir():
    raise FileNotFoundError(f"Folder not found: {model_path}")

# Process all XML files in folder and subfolders
process_xml_folder(model_path)
