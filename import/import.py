import logging
import sys

from institutions import prague_2lf_motol

DESTINATION_FOLDER = './data'


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    institution = sys.argv[1] if len(sys.argv) > 1 else input('Institution: ')
    input_folder = sys.argv[2] if len(sys.argv) > 2 else input('Path to the folder to import: ')

    if institution == 'prague-2lf-motol':
        prague_2lf_motol.import_patients(input_folder, DESTINATION_FOLDER)
    else:
        raise KeyError(f'The institution "{institution}" is unknown.')
