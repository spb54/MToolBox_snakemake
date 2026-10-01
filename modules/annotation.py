#!/usr/bin/env python
""" Functional annotation of VCFs with mtoolnote.

mtoolnote opens its SQLite database (inside the installed package) with the
default locking, which can hang on network filesystems where file locks are
unreliable. The database is only ever read, so open it read-only and
immutable, which needs no locks.
"""
import sys

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker


def use_lock_free_db():
    """ Make mtoolnote's human annotator open its database without locking. """
    from mtoolnote import human
    db = human.HumanAnnotator._dbfile
    engine = create_engine("sqlite:///file:{}?mode=ro&immutable=1&uri=true".format(db))
    human.HumanAnnotator.engine = engine
    human.HumanAnnotator.Session = sessionmaker(bind=engine)


def annotate_vcf(input_vcf, output_vcf, species="human"):
    """ Annotate input_vcf with mtoolnote, writing output_vcf. """
    import mtoolnote
    if species == "human":
        use_lock_free_db()
    mtoolnote.annotate(input_vcf, output_vcf, species=species)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("Usage: python -m modules.annotation input.vcf output.vcf species")
    annotate_vcf(*sys.argv[1:])
