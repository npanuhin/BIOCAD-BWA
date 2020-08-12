for /D /r %%a in (*) DO (
    if NOT "%%a" == "small/" (
        cd "%%a"
        cd "history/"
        ffmpeg -v warning -f image2 -i "%%d.png" -vf "fps=10,scale=1000:-1:flags=lanczos,palettegen" -y "palette.png"
        ffmpeg -v warning -framerate 2 -i "%%d.png" -i "palette.png" -lavfi "scale=1000:-1:flags=lanczos [x]; [x][1:v] paletteuse" -y "../history.gif"
        del "palette.png"
    )
)