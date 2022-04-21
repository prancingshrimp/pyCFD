#!/usr/bin/env python3.3

# from wand.image import Image
# from wand.display import display
# from wand.drawing import Drawing
import os
import shutil
from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw


path_to_files = "./pics_for_presentation"
resultsDir = path_to_files + "/" + "results"

templateImagePath = path_to_files + "/" + "animation_template.png"
templateImage = Image.open(templateImagePath)

field_list = ["u0", "v0", "p0"]
length = 964
width = 166
# ft = ImageFont.load("timR24.pil")

iter_step_list = []
for files in os.listdir(path_to_files + "/" + field_list[0]):
    if files.endswith(".png"):
        files = files[3:]
        files = files[:-4]
        if int(files) not in iter_step_list:
            iter_step_list.append(int(files))
iter_step_list.sort()

iter_step_list = [0]

# if os.path.exists(resultsDir) == True:
#     shutil.rmtree(resultsDir)
# os.mkdir(resultsDir)
for item in iter_step_list:
    print(item)
    for field in field_list:
        currentFile = path_to_files + "/" + field + "/" + field + "_" + str(item) + ".png"
        curImage = Image.open(currentFile)
        cutX = 30
        cutY = 301
        cutBox = (cutX, cutY, cutX+length, cutY+width)
        region = curImage.crop(cutBox)

        if field == field_list[0]:
            x = 0
            y = 54
        if field == field_list[1]:
            x = 0
            y = 291
        if field == field_list[2]:
            x = 0
            y = 530
        pasteBox = (x, y, x+length, y+width)
        templateImage.paste(region, pasteBox)
    savePath = path_to_files + "/" + "results" + "/" + "ani" + str(item) + ".png"
    print(savePath)
    templateImage.save(savePath)
    number = item * 10
    label = str(number)
    if number < 1000:
        label = "0" + str(number)
    if number < 100:
        label = "00" + str(number)
    if number == 0:
        label = "0000"
    os.system("convert " + savePath + " -font Helvetica -pointsize 30 -annotate +140+757 " + label + " " + savePath[:-4] + "_.png")
    os.remove(savePath)
    os.rename(savePath[:-4] + "_.png", savePath)
