import os
import subprocess
from performance import measure

input_folder = "automate_work"

@measure
def run_pixel_grower(file_path, angle, mode, inverse=False):
    subprocess.run(["./pixel_grower", file_path, str(1988), str(1988), str(400) ,str(angle), mode, "" if not inverse else "--inverse"])

for root, dirs, files in os.walk(input_folder):
    for file in files:
        file_path = os.path.join(root, file)
        for angle in [1, 2, 3, 5, 8, 10, 50, 100, 500, 1000, 5000, 10000]:
            for mode in ["avg", "min", "max"]:
                print(f"---- Running pixel grower for {file} with angles: {angle} and mode: {mode}. ----")
                run_pixel_grower(file_path, angle, mode, False)
                print("---- Normal done. ----")
                run_pixel_grower(file_path, angle, mode, True)
                print("---- Inverse done. ----")