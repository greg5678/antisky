# antisky
A github mirror of antisky

Videocrypt Cryptanalysis
------------------------

Markus Kuhn, University of Erlangen -- 1994-07-11

The files in this directory are TV pictures that have been received
from the ASTRA 1A satellite, transponder 8, station "Sky One". These
pictures have all been taken from a Star Trek TNG episode. The
original files are:

  d-raw.ppm       A picture without encryption showing "Deana Troi"
  d-vc1.ppm       The same scene (~2 frames later) with encryption
  e-vc1.ppm       An encrypted starship in space
  r-raw.ppm       Commander Riker on the bridge, not encrypted
  r-vc1.ppm       Same picture (~2 frames earlier) encrypted
  r-vc2.ppm       Same picture (~4 frames earlier) encrypted

Some of these files are available here in JPEG compressed format. Use
the djpeg program from <http://www.ijg.org/> to decompress these into
PPM. For other file formats, try <http://netpbm.sourceforge.net/>.

The images contain a PAL half frame (lines doubled) in an RGB color
space and were converted with a frame grabber (Parallax XVideo Board)
on a Sun SPARCstation 10.

The grayscale file 

  r-vc1.decrypted.pgm

demonstrates the results of an Videocrypt prototype cryptoanalysis
algorithm which I have developed. The only input to this algorithm was
the file r-vc1.ppm when this result was produced, NO information from
SmartCards or Videocrypt clone chips/cards was used. The algorithm is
still under development. Suggestions for improved algorithms are
welcome!

The file antisky.c contains the ANSI C source of the program that I
used to decrypt the image. E.g.

  antisky -1 -r20 r-vc1.ppm r-vc1.decrypted.pgm

produced the example file (~7 seconds on a Sun SPARCstation 10). Try
also

  antisky -1 -r20 -bc r-vc1.ppm r-vc1.decrypted.pgm

in order to see how the image looks after the first algorithm. The
additional line you'll see (option -c) is the edge formed by the
left/right border of the image which is detected by the second
algorithm. The first algorithm (cross-correlation) matches two
consecutive lines and shifts the lower line so that it matches the
next higher one. The second algorithm searches the border so that the
drift of the first algorithm and the total offset (which is inherited
by the first line) can be corrected.

Some theory:

Videocrypt rotates individual lines, or in other words, every line is
cut at a secret point in two parts and then both parts are
exchanged. I.e. if an original line in the pixtures was

  0123456789

(each digit represents one pixel), then the rotated version (here with
offset 3) looks like

  7890123456

What the first step of the ANTISKY algorithm is doing is only to
compare this rotated line in all 10 offsets

  7890123456
  6789012345
  5678901234
  ...
  9012345678
  8901234567

with the previous line. The measure of how good this line compares in
one particular offset to the previous one is the sum of the products
between pixels in the same column. In the output picture, consecutive
lines are rotated relative to each other, so that this measure is
maximized. The first line is not touched.

The general idea is quite simple, the question is how to implement
this idea as efficient as possible. The technical details of the
implementation are a little bit tricky and you won't understand them
without some mathematical background (fast fourier transform, real
FFT, cross correlation, convolution -> simple multiplication in
freqency domain, zero padding, etc.). If you want to learn, how all
this works, then definitely read the chapters 12.0-12.3 and 13.0-13.2
in the book

  William H. Press, Saul A. Teukolsky, William T. Vetterling:
  Numerical Recipes in C : The Art of Scientific Computing.
  Second edition, Cambridge University Press, 1993,
  ISBN 0521431085.
  <http://www.nr.com/>

or the chapters about FFT, convolution, and cross-correlation in any
introductory book about digital signal processing. It is really worth
to invest some time into studying this. I believe the Fast Fourier
Transform (FFT) is one of the coolest and most useful algorithms ever
invented.

Some care has to be taken prior to this algorithm, because frames
received from a frame grabber often look like this:

  @@789012345678@        (@ is a black pixel)

i.e. they have additional black borders and the left and right ends of
the line might overlap. Before the ANTISKY algorithm can produce good
results, you have to cut off black and doubled pixels (e.g. 3 from the
left and 2 from the right with options -l3 -r2). If you don't remove
these additional pixels, the corrected line would look like

  012345678@@@789

i.e. you get short black lines in the middle of the decoded image.

After the first step, the image has 2 defects (which you can see with
option -c):

  - the overall offset is the (unknown) offset of the first line
  - there are sometimes between 0 and ~3 pixels error in the offset, so the
    whole picture slowly drifts to the left or right.

In order to compensate these defects, the second step of the ANTISKY
algorithm searches the edge formed by the left and right border of the
image which is not visible in a correctly aligned image. The lines are
then rotated such that this edge vanishes. The edge detector algorithm
uses a dynamic programming technique, that is it calculates the
cheapest path from the top line down to any pixel based on the
cheapest path to the pixels in the line above. A path is "cheaper"
here if it goes along a high-contrast edge and if is does not deviate
too much from a vertical line (see the code for the exact cost
function used, I do not claim this one is optimal in any sense). If
you have no idea, what "dynamic programming" graph search algorithms
are, then read first of all the relevant chapter in
Cormen/Leiserson/Rivest: "Introduction to Algorithms", MIT Press,
1991, ISBN 0-262-03141-8.

It is not easy to find the correct edge which represents the image
margins in a picture, because some pictures contain many edges running
nearly vertically over the screen and some of them have often better
contrast than the left/right margin. But fortunately, a special
property of the videocrypt system allows to find the correct edge
quite safely. The following idea is implemented since antisky 0.92. If
you use option -m in order to make the cutting points in the image
visible, you'll discover, that cutting points appear never nearer than
<mindist> pixels to the image border and <mindist> is about 12% of the
image width. So the edge detection algorithm has to avoid getting
nearer than <mindist> pixles to any cutting point which excludes many
of the "alternative edges" that have confused previous releases of the
algorithm. You can make these deadzones around cutting points visible
with option -d.

Option -f allows you to reduce the horizontal resolution of the
algorithm to 1/2. This doubles the speed, but reduces image
quality. You can apply -f several times, e.g. -fff reduces the
resolution to 1/8 (which causes already severe distortions). Another
method of increasing the speed dramatically is to make the source
image smaller. An image that is obtained from the source image by
taking only every 4th pixel in horizontal and vertical direction
(i.e. which has only 1/16 of the original number of pixels) boosts the
decoders speed so that real-time TV is possible on good workstation
clusters (e.g. our fastest HP system here needs only 0.2s for such a
frame, so with 10 workstations you have nearly 50 Hz :-).

TODO:

  - Faster cross-correlation. This consumes nearly all the computing
    time at the moment and the FFT approach is quite well optimized
    since antisky 0.9. Anyone with better ideas (hierarchical
    matching, DSP, FFT processors, ...?).

  - Better edge detection (I have some ideas, but have not yet tried
    them). A lot of experimentation has to be done here.

  - Write drivers for the real time video recording harddisk cluster
    of the local computer graphics department. It might be used to
    decode a whole film during a night on a workstation cluster and
    put it back on video tape.

  - Color. The current software already supports it, but this is
    useless with our frame grabber (and every frame grabber I've seen
    so far). Problem: PAL decoder average the color vectors of the
    lines and if consecutive lines haven't similar colors, you'll get
    a random gray. A modified frame grabber hardware that doesn't
    implement the PAL standard color correction would be necessary!
    Experimental code for pure software PAL color correction is in the
    software since version 0.9, but it didn't work and thus has been
    deactivated with an #if 0. If you want to play around with the
    idea of doing the entire PAL decoding in software, you'll find on
    <http://www.ucl.ac.uk/~ucapwas/vd/paldec.html> a suitable code
    written by William Andrew Steer <ucapwas@ucl.ac.uk>. [That web
    page vanished in 2002, but an archival copy can be found, for
    example, on <http://web.archive.org/>].

A few more words about the color problem with Videocrypt, PAL and your
frame grabber:

The RGB video signal from the camera source is separated in the PAL
system in a luminance and a color (chrominance) part. Unfortunately,
every standard PAL frame grabber is confused by the videocrypt line
rotation and it destroys the chrominance part. One tricky detail of
the PAL system is that the chrominance parts of two lines are averaged
in order to avoid the phase-color problems (green faces, etc.) as they
appear in NTSC. But in Videocrypt, pairs of lines don't fit together,
so the averaged colors in the RGB image you get from your frame
grabber are always very chaotic. In order to allow antisky to compare
the lines correctly, only the part of the video signal which has
survived the frame grabber's PAL decoder, i.e. the luminance signal,
may be used. The luminance signal is what you see on an old
black/white TV set. According to the CCIR PAL standard, the luminance
Y of an RGB pixel is

    Y = R * 0.299 + G * 0.587 + B * 0.114.

The grayscale image that you get with these Y values should now be
undisturbed if your frame grabber strictly conforms to the PAL
standard. If your frame grabber allows you to store the image in YUV
format, then you should put only the Y component in a PGM file for
antisky as the frame grabber gives you already separated the luminance
Y and chrominance (U,V). If you give PPM RGB files to antisky, the
program will automatically extract the Y signal using the above
formula.

Markus

http://www.cl.cam.ac.uk/~mgk25/tv-crypt/image-processing/
