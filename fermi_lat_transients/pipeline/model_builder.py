"""
XML source model building for Fermi-LAT analysis.

Provides functions to create, modify, and manage XML source models
for likelihood analysis.
"""

import xml.etree.ElementTree as ET
from pathlib import Path


def build_source_model(catalog_fits, ra, dec, roi, galdiff, isodiff,
                       output_xml, make4FGLxml_path=None,
                       src_name=None, extended_dir=None):
    """
    Build an XML source model from the 4FGL catalog using make4FGLxml.

    Parameters
    ----------
    catalog_fits : str
        Path to 4FGL catalog FITS file.
    ra, dec : float
        ROI center coordinates [degrees].
    roi : float
        ROI radius [degrees].
    galdiff : str
        Path to Galactic diffuse model FITS.
    isodiff : str
        Path to isotropic diffuse template.
    output_xml : str
        Output XML model file path.
    make4FGLxml_path : str, optional
        Path to make4FGLxml.py script. If None, tries to import it.
    src_name : str, optional
        Name of any source to exclude from the model.
    extended_dir : str, optional
        Directory containing extended source templates.

    Returns
    -------
    str
        Path to the output XML model file.
    """
    import subprocess
    import sys

    if make4FGLxml_path is not None:
        # Run as external script
        cmd = [
            sys.executable, str(make4FGLxml_path),
            str(catalog_fits),
            str(output_xml),
            f"--galfile={galdiff}",
            f"--isofile={isodiff}",
            f"--ra={ra}", f"--dec={dec}",
            f"--radius={roi}",
        ]
        if extended_dir:
            cmd.append(f"--extDir={extended_dir}")
        subprocess.run(cmd, check=True)
    else:
        # Try to use make4FGLxml as a module
        try:
            from make4FGLxml import CatalogToXml
            catalog = CatalogToXml(
                catalog_fits, output_xml,
                galfile=galdiff, isofile=isodiff,
                ra=ra, dec=dec, radius=roi,
            )
            catalog.run()
        except ImportError:
            raise ImportError(
                "make4FGLxml not found. Provide make4FGLxml_path parameter "
                "or ensure make4FGLxml.py is importable."
            )

    return str(output_xml)


def add_transient_source(xml_path, name, ra, dec, spectrum='PowerLaw2',
                         output_xml=None, emin=100.0, emax=100000.0,
                         index=-2.0, flux=1.0e-7):
    """
    Add a transient point source to an existing XML model.

    Parameters
    ----------
    xml_path : str
        Path to existing XML model file.
    name : str
        Source name (e.g., 'Transient').
    ra, dec : float
        Source position [degrees].
    spectrum : str
        Spectral model type ('PowerLaw2' or 'PowerLaw').
    output_xml : str, optional
        Output file path. If None, modifies xml_path in-place.
    emin, emax : float
        Energy bounds [MeV] (for PowerLaw2).
    index : float
        Photon index (negative).
    flux : float
        Initial integral flux [ph/cm2/s].

    Returns
    -------
    str
        Path to the output XML file.
    """
    tree = ET.parse(xml_path)
    root = tree.getroot()

    # Build source element
    source = ET.SubElement(root, 'source', {
        'name': name,
        'type': 'PointSource',
    })

    if spectrum == 'PowerLaw2':
        spec = ET.SubElement(source, 'spectrum', {'type': 'PowerLaw2'})
        # Integral flux
        ET.SubElement(spec, 'parameter', {
            'name': 'Integral',
            'scale': '1.0',
            'value': str(flux),
            'min': '0.0',
            'max': '1000.0',
            'free': '1',
        })
        # Index
        ET.SubElement(spec, 'parameter', {
            'name': 'Index',
            'scale': '1.0',
            'value': str(index),
            'min': '-5.0',
            'max': '0.0',
            'free': '0',
        })
        # LowerLimit
        ET.SubElement(spec, 'parameter', {
            'name': 'LowerLimit',
            'scale': '1.0',
            'value': str(emin),
            'min': '20.0',
            'max': '1000000.0',
            'free': '0',
        })
        # UpperLimit
        ET.SubElement(spec, 'parameter', {
            'name': 'UpperLimit',
            'scale': '1.0',
            'value': str(emax),
            'min': '20.0',
            'max': '1000000.0',
            'free': '0',
        })
    else:
        # Simple PowerLaw
        spec = ET.SubElement(source, 'spectrum', {'type': 'PowerLaw'})
        ET.SubElement(spec, 'parameter', {
            'name': 'Prefactor',
            'scale': '1e-9',
            'value': str(flux * 1e9),
            'min': '0.0',
            'max': '1000.0',
            'free': '1',
        })
        ET.SubElement(spec, 'parameter', {
            'name': 'Index',
            'scale': '-1.0',
            'value': str(abs(index)),
            'min': '0.0',
            'max': '5.0',
            'free': '0',
        })
        ET.SubElement(spec, 'parameter', {
            'name': 'Scale',
            'scale': '1.0',
            'value': '100.0',
            'min': '50.0',
            'max': '500000.0',
            'free': '0',
        })

    # Spatial model (point source)
    spatial = ET.SubElement(source, 'spatialModel', {'type': 'SkyDirFunction'})
    ET.SubElement(spatial, 'parameter', {
        'name': 'RA', 'value': str(ra), 'free': '0',
        'scale': '1.0', 'min': '-360.0', 'max': '360.0',
    })
    ET.SubElement(spatial, 'parameter', {
        'name': 'DEC', 'value': str(dec), 'free': '0',
        'scale': '1.0', 'min': '-90.0', 'max': '90.0',
    })

    if output_xml is None:
        output_xml = xml_path
    tree.write(str(output_xml))
    return str(output_xml)


def remove_source_from_xml(xml_path, source_name, output_xml=None):
    """
    Remove a source from an XML model file.

    Parameters
    ----------
    xml_path : str
        Path to XML model file.
    source_name : str
        Name of the source to remove.
    output_xml : str, optional
        Output path. If None, modifies in-place.

    Returns
    -------
    str
        Path to the output XML file.
    """
    tree = ET.parse(xml_path)
    root = tree.getroot()

    for source in root.findall('source'):
        if source.get('name') == source_name:
            root.remove(source)
            break

    if output_xml is None:
        output_xml = xml_path
    tree.write(str(output_xml))
    return str(output_xml)


def freeze_all_in_xml(xml_path, output_xml=None, except_sources=None):
    """
    Set all parameters in an XML model to free='0' (frozen),
    optionally excluding specific sources.

    Parameters
    ----------
    xml_path : str
        Path to XML model file.
    output_xml : str, optional
        Output path. If None, modifies in-place.
    except_sources : list of str, optional
        Source names to keep free.

    Returns
    -------
    str
        Path to the output XML file.
    """
    if except_sources is None:
        except_sources = []

    tree = ET.parse(xml_path)
    root = tree.getroot()

    for source in root.findall('source'):
        if source.get('name') in except_sources:
            continue
        for param in source.iter('parameter'):
            param.set('free', '0')

    if output_xml is None:
        output_xml = xml_path
    tree.write(str(output_xml))
    return str(output_xml)
