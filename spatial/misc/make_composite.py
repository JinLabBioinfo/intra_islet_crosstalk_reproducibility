# ©2024 Bruker Spatial Biology, Inc. All rights reserved. Subject to additional license terms and conditions provided separately by Bruker Spatial Biology, Inc.
# https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/composite-images/making-composite-images.html
# Adapted from https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/blob/Main/assets/make-composite/make_composite.py

from PIL import Image, ImageSequence, ImageOps, ImageChops
import glob
import os
import re
import argparse
import numpy as np


def force_8bit(image):
    """Normalize image to 0–255 and convert to 8-bit grayscale (PIL 'L')."""
    array = np.array(image)
    array = array - array.min()
    if array.max() == 0:
        reduced_bit = 0
    else:
        reduced_bit = (array / array.max()) * 255
    return Image.fromarray(reduced_bit.astype("int8")).convert("L")


def compose_image(image_files, colors, channel_list, blackpoint, whitepoint):
    """Compose selected channels with screen blending; last file is used as the base layer."""
    img = ImageOps.colorize(
        Image.open(image_files[-1]),
        black="black",
        white=colors[-1],
        blackpoint=blackpoint,
        whitepoint=whitepoint,
        midpoint=np.floor((whitepoint + blackpoint * 2) / 3),
    )
    for image_file, color in zip(image_files[::-1][1:], colors[::-1][1:]):
        if image_file.split("_")[1] not in channel_list:
            continue
        img = ImageChops.screen(
            img,
            ImageOps.colorize(
                Image.open(image_file),
                black="black",
                white=color,
                blackpoint=blackpoint,
                whitepoint=whitepoint,
                midpoint=np.floor((whitepoint + blackpoint * 2) / 3),
            ),
        )
    return img


if __name__ == "__main__":
    # CLI: black/white points for contrast and channel filter (e.g., '01234' -> ch0..ch4)
    parser = argparse.ArgumentParser(...)
    parser.add_argument("--blackpoint", type=int, default=5)
    parser.add_argument("--whitepoint", type=int, default=250)
    parser.add_argument("--channel", type=str, default="01234")

    # Prepare output folders (expects empty/nonexistent dest dirs)
    morphology_dir = os.getcwd()
    os.chdir(morphology_dir)
    os.mkdir("raw_8bit")
    os.mkdir("composite")

    # Pass 1: explode multi-channel TIFFs into 8-bit per-channel images under raw_8bit/
    for image_file in glob.glob("*.tif") + glob.glob("*.TIF"):
        if re.compile(
            r"([\w]+_[\w]+_[\w]+_[\w]+_[\w]+_[\w]+_[F]+[\d]*)", flags=re.IGNORECASE
        ).match(image_file):
            fov_num = image_file.split("_")[-1].replace(".TIF", "")
            os.chdir(morphology_dir)
            image = Image.open(image_file)
            os.chdir(os.path.join(morphology_dir, "raw_8bit"))
            for channel, layer in enumerate(ImageSequence.Iterator(image)):
                force_8bit(layer).save(f"{fov_num}_ch{channel}_raw.tif")

    # Pass 2: compose selected channels per FOV into JPEG under composite/
    for fov in list(set([fov.split("_")[0] for fov in glob.glob("*.tif")])):
        os.chdir(os.path.join(morphology_dir, "raw_8bit"))
        image = compose_image(
            [image_file for image_file in glob.glob(fov + "*")],
            ["cyan", "red", "yellow", "blue", "magenta"],
            ["ch" + ch for ch in list(parser.parse_args().channel)],
            parser.parse_args().blackpoint,
            parser.parse_args().whitepoint,
        )
        os.chdir(os.path.join(morphology_dir, "composite"))
        image.save(f"CellComposite_{fov}.jpg", quality=95)
