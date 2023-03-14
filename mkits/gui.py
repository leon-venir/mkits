#from ase_notebook import AseView, ViewConfig, get_example_atoms




"""
from IPython.display import HTML
import ase.build
import ase.io


def atoms_to_html(atoms):
    'Return the html representation the atoms object as string'

    from tempfile import NamedTemporaryFile

    with NamedTemporaryFile('r+', suffix='.html') as ntf:
        atoms.write(ntf.name, format='html')
        ntf.seek(0)
        html = ntf.read()
    return html


def gui_show_ase_crys(inp:str):
    inpcell = ase.io.read(inp)
    tbut = ase.build.bulk(inpcell) 
    tbut_html = atoms_to_html(tbut)
    HTML(tbut_html)

"""

