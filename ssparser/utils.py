import re
import os

# return the regex to extract the pdb and metadata pkl files for each output
def alphafold_outfile_regex(multimer=False, relaxed=True):
    if relaxed:
        return (
            r"relaxed_model_.*_multimer_.*_pred_.*.pdb"
            if multimer
            else r"relaxed_model_.*_pred_.*.pdb"
        )
    else:
        return (
            r"unrelaxed_model_.*_multimer_.*_pred_.*.pdb"
            if multimer
            else r"unrelaxed_model_.*_pred_.*.pdb"
        )

# returns the .pkl and .pdb output of an AlphaFold attempt
def get_full_runs(data_dir: str, multimer=False, relaxed=True) -> tuple:
    re_str = alphafold_outfile_regex(multimer=multimer, relaxed=relaxed)
    r = re.compile(re_str)
    pairs = []

    for model in list(filter(r.match, os.listdir(data_dir))):
        run_output = model.strip(".pdb").strip("relaxed").strip("unrelaxed")
        pairs.append(
            (
                model,
                [
                    f
                    for f in os.listdir(data_dir)
                    if run_output in f and f.endswith(".pkl")
                ][0],
            )
        )
    return pairs


# get the .pdb files of each output
def get_run_pdb(data_dir: str, multimer=False, relaxed=True) -> list:
    re_str = alphafold_outfile_regex(multimer=multimer, relaxed=relaxed)
    return list(filter(re.compile(re_str).match, os.listdir(data_dir)))

# get all pdb files in a data directory
def get_all_pdb(data_dir: str) -> list:
    return [f for f in os.listdir(data_dir) if f.endswith(".pdb")]

# get all pdb files in a data directory
def get_all_csv(data_dir: str) -> list:
    return [f for f in os.listdir(data_dir) if f.endswith(".csv")]

def t_or_f(arg):
    ua = str(arg).upper()
    if "TRUE".startswith(ua):
        return True
    elif "FALSE".startswith(ua):
        return False
    else:
        raise BaseException(f"Argument {arg} is not a boolean value")
