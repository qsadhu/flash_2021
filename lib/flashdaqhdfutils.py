"""
    flashdaqhdfutils

    author:: Christopher Passow, Erland Muller
    copyright:: Christopher Passow 2020
"""
import glob
import sys
from datetime import datetime
from functools import reduce
from pathlib import Path
from typing import List, Iterable, Optional, Dict, Union

import h5py
from pandas import Series, DataFrame


assert sys.version_info >= (3, 8), f"Requires at least Python version 3.8, but was {sys.version}"


def dset_data_frame(group: h5py.Group, alias: Optional[str] = None) -> DataFrame:
    """Returns a pandas DataFrame for the given HDF group"""
    if not alias:
        alias = "/".join(group.name.split("/")[-3:])
    train_id = Series(group["index"], name="Train ID")
    np_array = group["value"][()] # decompress only once
    return Series((np_array[i] for i in train_id.index), name=alias, index=train_id).to_frame()


def all_value_dset_group_names(filename_pattern: str) -> List[str]:
    """ Returns a list of all value datasets in the files matching filename_pattern"""
    result = set()

    def value_dset_visitor(name, node: h5py.Group) -> None:
        if isinstance(node, h5py.Dataset) and name.endswith("/value"):
            result.add(node.parent.channel_name)

    for filename in glob.glob(filename_pattern):
        with h5py.File(filename, 'r') as h5:
            h5.visititems(value_dset_visitor)

    return sorted(result)


def run_files(run: int, *, directory: Union[Path, str], daq: str) -> List[Path]:
    """ Returns all files of given run located in directory for the given daq. """
    stream_name_prefixes = {"pbd": "GMD_DATA_gmd_data",
                            "pbd2": "FL2PhotDiag_pbd2_gmd_data",
                            "fl1user1": "FLASH1_USER1_stream_2",
                            "fl1user2": "FLASH1_USER2_stream_2",
                            "fl1user3": "FLASH1_USER3_stream_2",
                            "fl2user1": "FLASH2_USER1_stream_2",
                            "fl2user2": "FLASH2_USER2_stream_2"}
    date_time_section = lambda filename: str(filename).split("_")[-1]
    return sorted(Path(directory).glob(f"{stream_name_prefixes[daq]}_run{run}_*.h5"), key=date_time_section)


def files_dsets_data_frame(files: Iterable[Union[Path, str]],
                           dset_names: Union[Iterable[str], Dict[str, Optional[str]]]) -> DataFrame:
    """ Returns concatenated DataFrame for multiple files
        :parameter files: a non empty iterator of HDF files to search for datasets
        :parameter dset_names: an iterable of dataset names or a dictionary of dataset names and aliases
    """
    data_frames = [_single_file_dsets_data_frame(each_file, dset_names) for each_file in files]
    assert data_frames, "Assertion: at least one file in files, but files was empty"
    return reduce(DataFrame.combine_first, data_frames)


def date_time_data_frame(files: Iterable[Union[Path, str]]) -> DataFrame:
    """ Returns the date/time data frame for the given files """
    data_frames = [_single_file_date_time_data_frame(each_file) for each_file in files]
    return reduce(DataFrame.combine_first, data_frames)


def _single_file_dsets_data_frame(file_path: Union[Path, str],
                                 dset_names: Union[Iterable[str], Dict[str, Optional[str]]]) -> DataFrame:
    """ Returns a pandas DataFrame constructed for the given file.
        The DataFrame contains the datasets from the iterable or dict in the order specified by dset_names
    """
    if not isinstance(dset_names, dict):
        dset_names = {each: None for each in dset_names}  # since Python 3.7 preserves order
    with h5py.File(file_path, 'r') as h5_file:
        valid_names = (each_name for each_name in dset_names if each_name in h5_file)  # filter
        data_frames = (dset_data_frame(h5_file[each], dset_names[each]) for each in valid_names)  # map
        return reduce(DataFrame.combine_first, data_frames, DataFrame())


def _single_file_date_time_data_frame(file: Union[Path, str]) -> DataFrame:
    with h5py.File(file, 'r') as h5file:
        date_time_group = h5file[f"/FL2/Timing/Bunch pattern/train index 2"]
        train_id = Series(date_time_group["index"], name="Train ID")
        return DataFrame(date_time_group["time"],
                         columns=["Date/Time"],
                         index=train_id).applymap(datetime.fromtimestamp)
