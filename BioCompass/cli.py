# -*- coding: utf-8 -*-

import click

from .BioCompass import download_hits
from .BioCompass import download_mibig
from .BioCompass import createdb

@click.group()
def main():
    pass

@main.command(name="download-hits")
@click.option('--outputdir', default='./', type=click.Path(exists=True),
        help="Path to save the NCBI clusters.")
@click.argument('mgbfile', type=click.Path(exists=True), nargs=-1)
        #help="Multigeneblast file containing NCBI references to be downloaded.")
def downloadHits(mgbfile, outputdir):
    """Download NCBI clusters listed t in multigeneblast file."""
    for f in mgbfile:
        download_hits(f, outputdir)

@main.command(name="download-MIBiG")
@click.option('--outputdir', default='./', type=click.Path(exists=True),
        help="Path to save the MIBig genbank files.")
@click.option('--version', type=unicode, default='1.3',
        help="Version of MIBiG to download.")
def downloadMIBiG(outputdir, version):
    """Download MIBiG gbk database."""
    download_mibig(outputdir, version=version)

@main.command(name="createdb")
@click.option('--gbkdb', type=click.Path(exists=True),
        help='Path to downloaded gbk files')
@click.option('--outputdb', default=None, type=unicode,
        help='Path to database to be created')
@click.option('--multigeneblastdir', default=None, type=click.Path(exists=True),
        help='Path to database to be created')
@click.argument('clusterlist', type=unicode, nargs=-1)
def cli_createdb(gbkdb, outputdb, multigeneblastdir, clusterlist):
    createdb(gbkdb, outputdb, multigeneblastdir, clusterlist)


if __name__ == "__main__":
    main()
