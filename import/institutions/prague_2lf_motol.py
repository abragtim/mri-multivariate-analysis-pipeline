import gzip
import os
import re
import shutil
from logging import getLogger
from pathlib import Path

logger = getLogger(__name__)


def import_patients(input_folder, target_folder):
    logger.info('Importing patients from Prague 2LF Motol data...')

    filenames_map = {
        'rT1.nii': 'T1.nii',
        'rdmean.nii': 'dmean.nii',
        'rfa.nii': 'fa.nii',
        'rkmean.nii': 'kmean.nii',
        'kfa.nii': 'kfa.nii',
        'rNODDI__ficvf.nii': 'NODDI_ficvf.nii',
        'rNODDI__fiso.nii': 'NODDI_fiso.nii',
        'rNODDI__fmin.nii': 'NODDI_fmin.nii',
        'rNODDI__kappa.nii': 'NODDI_kappa.nii',
        'rNODDI__odi.nii': 'NODDI_odi.nii',
        'c1rT1.nii': 'GMSegmentation.nii',
        'c2rT1.nii': 'WMSegmentation.nii',
        'Segmentation-T1W-label.nii': 'lesionSegmentation.nii',
        'Segmentation-T1W-label_1.nii': 'lesionSegmentation.nii',
        'Segmentation-rT1W-label.nii': 'lesionSegmentation.nii',
    }
    regex_filenames_keys = [fr'.*{key}.*' for key in filenames_map
                            if key != 'rT1.nii']

    def map_filename(input_filename):
        target_filename = filenames_map.get(input_filename)
        if target_filename:
            return target_filename

        matched_keys = [regex_key.replace('.*', '')  # extract pattern back to find in filenames_map
                        for regex_key in regex_filenames_keys
                        if re.match(regex_key, input_filename)]
        if not matched_keys:
            return None
        elif len(matched_keys) == 1:
            return filenames_map[matched_keys[0]]
        else:
            raise IndexError(f"More than one regex key matched for the patient: {matched_keys}. "
                             "This is internal error, fix the code.")

    patients = os.listdir(input_folder)
    for patient in patients:
        logger.info(f'Importing patient {patient}...')
        patient_input_path, patient_target_path = _get_patient_paths(patient, input_folder,
                                                                     target_folder)

        expected_files = set(filenames_map.values())
        for root, _, files in os.walk(patient_input_path):
            for file in files:
                file_input_path = Path(f'{root}/{file}')
                if file.endswith('.gz'):
                    file_input_path = _extract_and_delete_gzip(str(file_input_path))
                    file = file.replace('.gz', '')

                target_file = map_filename(file)
                if not target_file:
                    continue
                file_target_path = patient_target_path / target_file

                shutil.copy(file_input_path, file_target_path)
                expected_files.discard(file_target_path.name)

        if expected_files:
            raise ValueError(f'Some files for the patient {patient} were not imported: '
                             f'{", ".join(expected_files)}.')

    logger.info('Successfully imported patients from Prague 2LF Motol data.')


def _get_patient_paths(patient, input_folder, target_folder):
    patient_input_path = Path(os.path.join(input_folder, patient))
    patient_target_path = Path(os.path.join(target_folder, patient))
    if not os.path.exists(patient_target_path):
        #os.mkdir(patient_target_path) # to create new folders, do not forget to add metadata.json manually
        raise ValueError(f'Target path {patient_target_path} does not exist. '
                         f'Please ensure that the patient folder is created within '
                         f'the data directory before proceeding with the import.')
    return patient_input_path, patient_target_path


def _extract_and_delete_gzip(filepath: str) -> str:
    extracted_file_path = filepath.replace('.gz', '')
    with gzip.open(filepath, 'rb') as gz_file:
        with open(extracted_file_path, 'wb') as extracted_file:
            logger.debug(f'Extracting {filepath} file...')
            shutil.copyfileobj(gz_file, extracted_file)

    logger.debug(f"Removing {filepath} file...")
    os.remove(filepath)

    return extracted_file_path
