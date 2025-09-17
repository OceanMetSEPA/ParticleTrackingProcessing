from pathlib import Path
from .xml2dict import xml2dict

def process_xml_folder(model_path):
    """
    Process all XML files in a folder (including subfolders).
    
    Parameters
    ----------
    model_path : str or Path
        Folder containing XML files.
    """

    # Ensure Path object
    model_path = Path(model_path)

    # Convert to string for legacy string checks (UNC paths etc.)
    model_path_str = str(model_path)

    # Fix leading single backslash if needed (for network paths)
    if model_path_str.startswith("\\") and not model_path_str.startswith("\\\\"):
        model_path_str = "\\" + model_path_str
        model_path = Path(model_path_str)

    if not model_path.is_dir():
        raise FileNotFoundError(f"Folder not found: {model_path}")

    # Recursively find all XML files
    xml_files = list(model_path.rglob("*.xml"))
    print(f"{len(xml_files)} XML files found:")
    print(xml_files)

    for xml_file in xml_files:
        matfile_name = str(xml_file).replace(".xml", ".mat")
        matfile_name = (
            matfile_name.replace("pt3D", "")
                        .replace("5minUnComp", "_trackStruct")
                        .replace("__", "_")
                        .replace("_ECLH", "")
                        .replace("_WLLS", "")
                        .replace("_FOC", "")
                        .replace("_ES", "")
        )

        if not Path(matfile_name).is_file():
            print(f"Processing {xml_file}")
            xmlDict, matfile_name = xml2dict(xml_file, matfileName=matfile_name)
        else:
            print(f"{matfile_name} exists already!")

    print("Done!")
