#%%
import subprocess
import glob
import os
from multiprocessing import Pool
#%%
class ObsDataSync:
    def __init__(self, source_base="/lyman/data1/obsdata", dest_base="/home/hhchoi1022/data/7DT/obsdata"):
        """
        Initializes the ObsDataSync class.
        
        :param source_base: Base source directory containing telescope folders (7DT??)
        :param dest_base: Base destination directory where data will be synced
        """
        self.source_base = source_base
        self.dest_base = dest_base

    def _run_rsync(self, telescope_id, folder_keys):
        """
        Runs rsync for a given telescope directory (7DT??) and folder keys.
        
        :param telescope_id: The specific telescope folder (e.g., '7DT01', '7DT02')
        :param folder_keys: List of folder keys to sync
        """
        source_dir = os.path.join(self.source_base, telescope_id)
        dest_dir = os.path.join(self.dest_base, telescope_id)

        for folder_key in folder_keys:
            src = os.path.join(source_dir, folder_key)
            dest = os.path.join(dest_dir, folder_key)
            if not os.path.exists(dest):  # Ensure destination exists
                os.makedirs(dest)
            
            if os.path.exists(src):  # Ensure source exists
                cmd = ["rsync", "-av", "--progress", src + "/", dest + "/"]
                print(f"Starting rsync for {src} -> {dest}")
                subprocess.run(cmd)
            else:
                print(f"Skipping missing source directory: {src}")

    def sync_all_folders(self, folder_keys):
        """
        Finds all matching 7DT?? folders and runs rsync processes in parallel.
        
        :param folder_keys: List of folder keys to sync
        """
        if isinstance(folder_keys, str):
            folder_keys = [folder_keys]
        source_pattern = os.path.join(self.source_base, "7DT??")
        telescope_dirs = glob.glob(source_pattern)  # Get list of matching directories

        telescope_ids = [os.path.basename(os.path.normpath(t)) for t in telescope_dirs]  # Extract 7DT??

        # Run rsync processes in parallel
        with Pool(processes=len(telescope_ids)) as pool:
            pool.starmap(self._run_rsync, [(tid, folder_keys) for tid in telescope_ids])
#%%
# Example usage:
if __name__ == "__main__":
    folder_keys = ["2025-02-02_gain2750"]  # Add required folder keys
    sync_manager = ObsDataSync()
    sync_manager.sync_all_folders(folder_keys)

# %%
