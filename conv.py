import cv2
import glob

# encoder(for mp4)
fourcc = cv2.VideoWriter_fourcc('m', 'p', '4', 'v')
# output file name, encoder, fps, size(fit to image size)
video = cv2.VideoWriter('video.mp4',fourcc, 30.0, (640, 480))
l = glob.glob('img/*.png')
l.sort()
for i in l:
    # hoge0000.png, hoge0001.png,..., hoge0090.png
    img = cv2.imread(i)
    # add
    video.write(img)
    print(i)

video.release()
print('written')