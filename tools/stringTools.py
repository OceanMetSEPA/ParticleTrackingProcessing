import os

def siteNameFromString(path: str) -> str:
    """
    Given a string/filename, try to work out the site name.

    Parameters
    ----------
    path : str
        Full path or filename

    Returns
    -------
    str
        Extracted site name
    """
    if not path:
        return ""

    # Extract final part of path without extension
    base = os.path.splitext(os.path.basename(str(path)))[0]

    # Strings to remove
    str2remove = [
        'MIKE2025','ES','dfs0_','PFOW_SSM_5MinPart_','_surfaceConc','_trackStruct',
        'TrackStruct','trackStruct','pt3D','5minUnComp','_dfsu','WK12ToEndOfMayDecoup_',
        '_5minPart','_Decoupled_','ECLH1993','FOC2019','FOC_','WLLS1993','_ECLH',
        '_PFOW','_','ES2025','2025'
    ]

    for s in str2remove:
        base = base.replace(s, "")

    return base


def modelDomainFromString(path: str) -> str:
    """
    Given a string/filename, try to determine the model domain.

    Parameters
    ----------
    path : str
        String to inspect

    Returns
    -------
    str
        Model domain ('' if none found)

    Raises
    ------
    ValueError
        If multiple model domains are found
    """
    model_domains = ['FOC', 'WLLS', 'ECLH', 'WestCOMS', 'PFOW', 'EastSkye']

    found = [d for d in model_domains if d in str(path)]

    if len(found) == 0:
        return ''
    elif len(found) == 1:
        return found[0]
    else:
        raise ValueError(f"{path} has more than one model domain in it!")
