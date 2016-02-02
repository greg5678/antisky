/*
 *  ANTISKY -- Videocrypt cryptoanalysis algorithm
 *
 *  This program reconstructs TV images that have been encrypted
 *  by pay-TV broadcasting stations that use the Videocrypt encryption
 *  system. The algorithm uses only the statistical properties of normal
 *  TV pictures and requires NO secret information that is stored
 *  in the Videocrypt SmartCards or their illegal clones.
 *
 *  This software is only provided for scientific demonstration
 *  purposes as an example of an interesting signal processing
 *  algorithm. The author rejects any responsibility for pay-TV
 *  piracy abuses of this software.
 *
 *  Markus Kuhn <mskuhn@cip.informatik.uni-erlangen.de>
 *
 *  $Id: antisky.c,v 1.5 1994/02/20 15:14:01 mskuhn Exp $
 */


#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#ifndef SEEK_CUR
#define SEEK_CUR 1
#endif

#define USE_FFT

typedef struct {
  FILE *f;
  enum {PGM_raw, PPM_raw} type;
  int xmax, ymax;
} G_FILE;

typedef unsigned char RGB_line[3];

float  *twr[2], *twi[2];     /* base vectors for FFT */

char version[] = "ANTISKY 0.92 -- Videocrypt cryptoanalysis algorithm"
                 " -- Markus Kuhn\n";

char usage[] = "%s\n%s [options] [<encrypted file> | -] [<decrypted file>]\n"
"\nThe encrypted file is a single frame of a video sequence that has been\n"
"encrypted using the Videocrypt system. The output file will be a b/w\n"
"PGM file with the decrypted image. If no filenames are given, stdin and\n"
"stdout are used.\n\ncommand line options:\n\n"
"  -r<number-of-columns>    The number of pixel columns that should be\n"
"                           removed from the right margin of the source\n"
"                           image. It is very important that this parameter\n"
"                           or -l is adjusted correctly, because the source\n"
"                           image might contain a few pixels doubled on\n"
"                           both sides of the image which will after\n"
"                           decoding be inserted somewhere in the middle\n"
"                           of the image and cause matching problems.\n"
"  -l<number-of-columns>    Number of pixel columns to be removed from the\n"
"                           left margin before decoding starts. If your\n"
"                           frame grabber has cut off too many pixels at\n"
"                           the borders, try negative values here.\n"
"  -v                       Verbose. Print some additional information.\n"
"  -b                       Mark border found by edge detection algorithm.\n"
"                           Best applyied together with -c.\n"
"  -m                       Mark the positions of the cutting points in the\n"
"                           encoded image (only for demonstration).\n"
"  -d                       Mark the deadzones near cutting points where\n"
"                           no image border can be (only for demonstration).\n"
"  -c                       Do only maximize the interline cross-correlation\n"
"                           (algorithm 1) and do not search the image margin\n"
"                           and correct it. (Only for demonstration of\n"
"                           how the two steps of the process work.)\n"
"  -o                       Use the input format instead of PGM for output\n"
"                           (currently also supported: PPM).\n"
"  -f                       Reduce the horizontal resolution to 50%%. This\n"
"                           makes the program faster, but reduces quality.\n"
"  -0                       Process all lines (default).\n"
"  -1                       Process only odd lines and double them in\n"
"                           the output file.\n"
"  -2                       Process only even lines and double them.\n"
"\n";


void *checkedmalloc(size_t n) {
  void *p;
  
  if ((p = malloc(n)) == NULL) {
    fprintf(stderr, "Sorry, not enough memory available!\n");
    exit(1);
  }
  
  return p;
}


/* generate vector tables for FFT */

void init_fft(int nn)
{
  int n, mmax, m, istep;
  double wtemp, wpr, wpi, theta;
  int c, v, isign;
  
  n=nn << 1;

  c = 0;
  mmax = 2;
  while (n > mmax) {
    istep = mmax << 1;
    c += mmax / 2 + 1;
    mmax = istep;
  }

  for (v = 0; v <= 1; v++) {
    twr[v] = (float *) checkedmalloc(c * sizeof(float));
    twi[v] = (float *) checkedmalloc(c * sizeof(float));
    isign = v ? -1 : 1;
    c = 0;
    mmax = 2;
    while (n > mmax) {
      istep = mmax << 1;
      theta = isign * (6.28318530717959 / mmax);
      wtemp = sin(0.5 * theta);
      wpr = -2.0 * wtemp * wtemp;
      wpi = sin(theta);
      twr[v][c]=1.0;
      twi[v][c]=0.0;
      c++;
      for (m = 1; m < mmax; m += 2) {
	twr[v][c] = twr[v][c-1] * wpr - twi[v][c-1] * wpi + twr[v][c-1];
	twi[v][c] = twi[v][c-1] * wpr + twr[v][c-1] * wpi + twi[v][c-1];
	c++;
      }
      mmax = istep;
    }
  }

  return;
}


/* 1-D real FFT implementation */

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void fft(float data[], int nn, int isign)
{
  int n, mmax, m, j, istep, i;
  float wr, wi;
  float tempr, tempi;
  float *ptwr, *ptwi;
  
  n = nn << 1;
  j = 1;
  for (i = 1;i < n;i += 2) {
    if (j > i) {
      SWAP(data[j], data[i]);
      SWAP(data[j+1], data[i+1]);
    }
    m = n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  ptwr = twr[isign < 0];
  ptwi = twi[isign < 0];
  mmax = 2;
  while (n > mmax) {
    istep = mmax << 1;
    wr = *(ptwr++);
    wi = *(ptwi++);
    for (m = 1; m < mmax; m += 2) {
      for (i = m; i <= n; i += istep) {
	j = i + mmax;
	tempr = wr*data[j]   - wi * data[j+1];
	tempi = wr*data[j+1] + wi * data[j];
	data[j]   = data[i]   - tempr;
	data[j+1] = data[i+1] - tempi;
	data[i]   += tempr;
	data[i+1] += tempi;
      }
      wr = *(ptwr++);
      wi = *(ptwi++);
    }
    mmax = istep;
  }
}

#undef SWAP


void real_fft(float data[], int n, int isign)
{
  int i, i1, i2, i3, i4, np3;
  float c1=0.5, c2, h1r, h1i, h2r, h2i;
  double wr, wi, wpr, wpi, wtemp, theta;
  
  theta = 3.141592653589793 / (double) (n>>1);
  if (isign == 1) {
    c2 = -0.5;
    fft(data, n>>1, 1);
  } else {
    c2 = 0.5;
    theta = -theta;
  }
  wtemp = sin(0.5 * theta);
  wpr = -2.0 * wtemp * wtemp;
  wpi = sin(theta);
  wr = 1.0 + wpr;
  wi = wpi;
  np3 = n+3;
  for (i = 2;i <= (n>>2); i++) {
    i4 = 1 + (i3 = np3 - (i2 = 1 + (i1 = i + i - 1)));
    h1r =  c1 * (data[i1] + data[i3]);
    h1i =  c1 * (data[i2] - data[i4]);
    h2r = -c2 * (data[i2] + data[i4]);
    h2i =  c2 * (data[i1] - data[i3]);
    data[i1] =  h1r + wr * h2r - wi * h2i;
    data[i2] =  h1i + wr * h2i + wi * h2r;
    data[i3] =  h1r - wr * h2r + wi * h2i;
    data[i4] = -h1i + wr * h2i + wi * h2r;
    wr = (wtemp = wr) * wpr -    wi * wpi + wr;
    wi =          wi  * wpr + wtemp * wpi + wi;
  }
  if (isign == 1) {
    data[1] = (h1r = data[1]) + data[2];
    data[2] = h1r - data[2];
  } else {
    data[1] = c1 * ((h1r = data[1]) + data[2]);
    data[2] = c1 * (h1r - data[2]);
    fft(data, n>>1, -1);
  }
}


/* graphic file handling routines */

/* skip whitespace and comments in PGM/PPM headers */
void skip_white(FILE *f)
{
  int c;
  
  do {
    while (isspace(c = getc(f)));
    if (c == '#')
      while ((c = getc(f)) != '\n' && c != EOF);
    else {
      ungetc(c, f);
      return;
    }
  } while (c != EOF);
  
  return;
}

/* read header of input file */
int g_openr(G_FILE *gf)
{
  char magic[10];
  int vmax;
  
  if (gf->f == NULL) abort();
  magic[0] = fgetc(gf->f);
  magic[1] = fgetc(gf->f);
  if (magic[0] == 'P' && isdigit(magic[1])) {
    /* PGM or PPM */
    fgets(magic + 2, 10, gf->f);
    if (!strcmp(magic, "P5\n")) gf->type = PGM_raw;
    else if (!strcmp(magic, "P6\n")) gf->type = PPM_raw;
    else {
      fprintf(stderr, "Unsupported input file type!\n");
      return(1);
    }
    skip_white(gf->f);
    fscanf(gf->f, "%d", &gf->xmax);
    skip_white(gf->f);
    fscanf(gf->f, "%d", &gf->ymax);
    skip_white(gf->f);
    fscanf(gf->f, "%d", &vmax);
    getc(gf->f);
  } else {
    fprintf(stdin, "Unsupported input file type!\n");
    return(1);
  }
  
  return 0;
}

/* write header of output file */
void g_openw(G_FILE *gf)
{
  switch (gf->type) {
  case PGM_raw:
    fprintf(gf->f, "P5\n# %s%d %d\n%d\n", version,
	    gf->xmax, gf->ymax, 255);
    break;
  case PPM_raw:
    fprintf(gf->f, "P6\n# ANTISKY -- Markus Kuhn\n%d %d\n%d\n",
	    gf->xmax, gf->ymax, 255);
    break;
  default: 
    abort();
  }

  return;
}

/* read a single image line and store it as RGB values */
void g_read_rgb(RGB_line *line, G_FILE *gf)
{
  int i;
  
  switch (gf->type) {
  case PPM_raw:
    if (line)
      fread(line, 3, gf->xmax, gf->f);
    else
      for (i = 0; i < 3 * gf->xmax; i++)
	getc(gf->f);
    break;
  case PGM_raw:
    if (line)
      for (i = 0; i < gf->xmax; i++)
	line[i][0] = line[i][1] = line[i][2] = getc(gf->f);
    else
      for (i = 0; i < gf->xmax; i++)
	getc(gf->f);
    break;
  default:
    abort();
  }
  
  return;
}

/* write a single line of RGB values to a file */
void g_write_rgb(RGB_line *line, int offset, G_FILE *gf)
{
  int i, v;
  RGB_line *p;
  
  switch (gf->type) {
  case PPM_raw:
    fwrite(line + offset, 3, gf->xmax - offset, gf->f);
    fwrite(line, 3, offset, gf->f);
    break;
  case PGM_raw:
    p = line + offset;
    for (i = 0; i < gf->xmax; i++) {
      v = (3 * (int) (*p)[0] + 
	   6 * (int) (*p)[1] + 
	   1 * (int) (*p)[2]) / 10;
      if (v < 0) v = 0;
      if (v > 255) v = 255;
      putc(v, gf->f);
      if (p++ == (line + gf->xmax)) p = line;
    }
    break;
  default:
    abort();
  }
  
  return;
}


int main(int argc, char **argv)
{
  char *fnin = NULL, *fnout = NULL;
  int line, lines, linecount = 0, i, j, k, x, fftsize, lumres;
  RGB_line **rgb;
  float *lum[2];  /* luminance signal of this and previous line */
  float *f[2];    /* fourier transformed lum */
  float *cc;      /* cross correlation between this and previous line */
#ifndef USE_FFT
  float *p, sum;
#else
  float dum;
#endif
  float maxsum = 0;
  int last = 1, new = 0;
  G_FILE gin, gout;
  int *offset, *hist;
  signed char **dir;   /* direction vectors for dynamic programming alg. */
  int *val[2];         /* bonus values for dyn. prog. algorithm */
  int mindist;         /* this is deadzone in lumres pixels */
  /* parameters */
  int from_stdin = 0;
  int rmargin = 0, lmargin = 0;
  int halfframe = 0;  /* 1: only first, 2: only second, 0: both halfframes */ 
  int verbose = 0;
  int markbreak = 0;
  int mark_edge = 0;
  int corr_edge = 1;
  int max_edge_bias = 4;
  int orig_fmt = 0;
  int resolution = 1;
  int show_deadzone = 0;
  float deadzone = 0.12;   /* this multiplied with the image width is the
			    * minimum distance between the real image's
			    * left/right border and any cutting-point */

  /* read command line arguments */
  gin.f = stdin;
  gout.f = stdout;
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '/' && argv[i][1] == '?' && argv[i][2] == 0) {
      fprintf(stderr, usage, version, argv[0]);
      exit(1);
    } else if (argv[i][0] == '-')
      if (argv[i][1] == '\0') {
        if (from_stdin) {
          fprintf(stderr, usage, version, argv[0]);
          exit(1);
        }
        from_stdin = 1;
      } else
        for (j = 1; j < 2000 && argv[i][j] != 0; j++)
          switch (argv[i][j]) {
          case 'r':
          case 'R':
            if (argv[i][j+1] != '\0') {
              rmargin = atoi(argv[i] + j + 1);
              j = 9999;
            } else {
              fprintf(stderr, usage, version, argv[0]);
              exit(1);
            }
            break;
          case 'l':
          case 'L':
            if (argv[i][j+1] != '\0') {
              lmargin = atoi(argv[i] + j + 1);
              j = 9999;
            } else {
              fprintf(stderr, usage, version, argv[0]);
              exit(1);
            }
            break;
          case 'b':
          case 'B':
            mark_edge = 1;
            break;
          case 'm':
          case 'M':
            markbreak = 1;
            break;
          case 'd':
          case 'D':
            show_deadzone = 1;
            break;
          case 'c':
          case 'C':
            corr_edge = 0;
            break;
          case 'o':
          case 'O':
            orig_fmt = 1;
            break;
          case 'f':
          case 'F':
            resolution *= 2;
            break;
          case '0':
          case '1':
          case '2':
            halfframe = argv[i][j] - '0';
            break;
          case 'v':
          case 'V':
            fprintf(stderr, "%s\n", version);
            verbose = 1;
            break;
          default:
            fprintf(stderr, usage, version, argv[0]);
            exit(1);
          }
    else if (gin.f == stdin && !from_stdin) {
      fnin = argv[i];
      gin.f = fopen(fnin, "rb");
      if (gin.f == NULL) {
        fprintf(stderr, "Can't open input file '%s", fnin);
        perror("'");
        exit(1);
      }
    } else if (gout.f == stdout) {
      fnout = argv[i];
      gout.f = fopen(fnout, "wb");
      if (gout.f == NULL) {
        fprintf(stderr, "Can't open output file '%s", fnout);
        perror("'");
        exit(1);
      }
    } else {
      fprintf(stderr, usage, version, argv[0]);
      exit(1);
    }
  }
  
  if (g_openr(&gin)) {
    if (fnin == NULL)
      fprintf(stderr, "Can't read file from stdin!\n");
    else
      fprintf(stderr, "Can't read input file '%s'!\n", fnin);
    exit(1);
  }
  gout.xmax = gin.xmax - lmargin - rmargin;
  gout.ymax = gin.ymax;
  gout.type = orig_fmt ? gin.type : PGM_raw;
  lines = gin.ymax;
  if (halfframe) {
    lines = gin.ymax / 2;
    gout.ymax = 2 * lines;
  }
  g_openw(&gout);
  lumres = gout.xmax / resolution;
  fftsize = 2;
  for (i = lumres - 1; i > 0; i >>= 1) fftsize <<= 1;
  mindist = gout.xmax * deadzone / resolution; 

  /* allocate buffers */
  rgb = (RGB_line **) checkedmalloc(sizeof (RGB_line *) * lines);
  dir = (signed char **) checkedmalloc(sizeof (signed char *) * lines);
  for (i = 0; i < lines; i++) {
    rgb[i] = (RGB_line *) checkedmalloc(sizeof (RGB_line) * gin.xmax);
    dir[i] = (signed char *) checkedmalloc(sizeof (signed char) * lumres);
  }
  offset = (int *) checkedmalloc(sizeof (int) * lines);
  f[last] = (float *) checkedmalloc(sizeof (float) * fftsize);
  f[new] = (float *) checkedmalloc(sizeof (float) * fftsize);
  cc = (float *) checkedmalloc(sizeof (float) * fftsize);
  lum[last] = (float *) checkedmalloc(sizeof (float) * 2 * lumres);
  lum[new] = (float *) checkedmalloc(sizeof (float) * 2 * lumres);
  hist = (int *) checkedmalloc(sizeof (int) * lumres);
  val[last] = (int *) checkedmalloc(sizeof (int) * lumres);
  val[new] = (int *) checkedmalloc(sizeof (int) * lumres);

  for (x = 0; x < lumres; x++) {
    hist[x] = 0;
    val[last][x] = 0;
    val[new][x] = 0;
  }
  
  /* load source image into rgb */
  if (verbose)
    fprintf(stderr, "loading data ...\n");
  for (line = 0; line < lines; line++) { 
    if (halfframe == 2) g_read_rgb(NULL, &gin);
    g_read_rgb(rgb[line], &gin);
    if (halfframe == 1) g_read_rgb(NULL, &gin);
  }

#if 0
  /* correct PAL color distortion -- this didn't work */
  float y, *u, *v;  /* PAL color components */
  u = (float *) checkedmalloc(sizeof (float) * gin.xmax);
  v = (float *) checkedmalloc(sizeof (float) * gin.xmax);
  if (lines > 0)
    for (x = 0; x < gin.xmax; x++) {
      u[i] = (- 0.147 * (float) rgb[0][x][0] 
	      - 0.288 * (float) rgb[0][x][1] 
	      + 0.436 * (float) rgb[0][x][2]);
      v[i] = (+ 0.615 * (float) rgb[0][x][0]
	      - 0.514 * (float) rgb[0][x][1]
	      - 0.101 * (float) rgb[0][x][2]);
    }
  for (line = 1; line < lines; line++) { 
    for (x = 0; x < gin.xmax; x++) {
      y = (+ 0.299 * (float) rgb[line][x][0]
	   + 0.587 * (float) rgb[line][x][1]
	   + 0.114 * (float) rgb[line][x][2]);
      u[x] = 2 * (- 0.147 * (float) rgb[line][x][0] 
		  - 0.288 * (float) rgb[line][x][1] 
		  + 0.436 * (float) rgb[line][x][2]) - u[x];
      v[x] = 2 * (+ 0.615 * (float) rgb[line][x][0]
		  - 0.514 * (float) rgb[line][x][1]
		  - 0.101 * (float) rgb[line][x][2]) - v[x];
      /* red */
      k = 1.14 * v[x] + y + 0.49;
      if (k < 0) k = 0; else if (k > 255) k = 255;
      rgb[line][x][0] = k;
      k = - 0.58 * v[x] - 0.39 * u[x]+ y + 0.49;
      if (k < 0) k = 0; else if (k > 255) k = 255;
      rgb[line][x][1] = k;
      k = 2.04 * u[x] + y + 0.49;
      if (k < 0) k = 0; else if (k > 255) k = 255;
      rgb[line][x][2] = k;
    }
  }
#endif

  /* main decryption loop */
  offset[0] = 0;
  init_fft(fftsize);
  for (line = 0; line < lines; line++) {
    for (i = 0; i < lumres; i++) {
      lum[new][i] = 0;
      for (j = 0; j < resolution; j++) {
	x = i * resolution + j + lmargin;
	if (x >= 0)
	  lum[new][i] +=
	    (0.299 * (float) rgb[line][x][0] +
	     0.587 * (float) rgb[line][x][1] +
	     0.114 * (float) rgb[line][x][2]);
      }
      lum[new][i] /= resolution;
    }
#ifdef USE_FFT
    /* transform luminance of each line into the frequency domain
     * using the fast fourier transform                            */
    for (i = lumres; i < fftsize; i++)
      f[new][i] = 0.0;   /* zero padding */
    memcpy(f[new], lum[new], sizeof(float) * lumres);
    real_fft(f[new] - 1, fftsize, 1);  /* FFT */
#endif
    if (line) {
#ifdef USE_FFT
      /* compute cross-correlation using real FFT convolution */
      /* complex multiplication in the frequency domain is
	 equivalent to convolution in time domain */
      f[last][0] *= f[new][0];    /* DC part and sample-frequency/2 */
      f[last][1] *= f[new][1];    /* are always real */
      for (i = 2; i < fftsize; i += 2) {
	/* multiplication of complex numbers */
	dum = f[last][i];
	f[last][i  ] = f[last][i] * f[new][i]    + f[last][i+1] * f[new][i+1];
	f[last][i+1] = dum        * f[new][i+1]  - f[last][i+1] * f[new][i];
      }
      real_fft(f[last] - 1, fftsize, -1);  /* inverse transform */
      /* find peak in cross-correlation */
      maxsum = -1e20;
      for (x = 0; x < lumres; x++)
	if (f[last][x] + f[last][fftsize - lumres + x] > maxsum) { 
	  maxsum = f[last][x] + f[last][fftsize - lumres + x];
	  offset[line] = x;
	}
#else
      /*
       * Simple two-loop cross-correlation. Mathematically equivalent to
       * above FFT code, but slower. This is only left in here for
       * testing and in order to help explaining the algorithm to
       * people who are not familiar with fourier transform, convolution, etc.
       */
      maxsum = -1e20;
      for (x = 0; x < lumres; x++) {
	sum = 0.0;
	p = lum[new] + x;
	for (i = 0; i < lumres - x; i++)
	  sum += lum[last][i] * *(p++);
	p = lum[new];
	for (; i < lumres; i++)
	  sum += lum[last][i] * *(p++);
	if (sum > maxsum) maxsum = sum, offset[line] = x;
      }
#endif
      offset[line] = (offset[line] + offset[line - 1]) % lumres;
      if (corr_edge || mark_edge) {
	/* dynamic programming edge search, part 1 */
	for (i = 0; i < lumres; i++) {
	  val[new][i] = 0;
	  j = i - max_edge_bias;
	  if (j < 0) j += lumres;
	  for (x = -max_edge_bias; x <= max_edge_bias; x++) {
	    k = x;             /* k is the penalty for drifting */
	    if (k < 0) k = -k;
	    if (k == 1) k = 0;
	    k *= 2 * resolution;
	    if (val[last][j]-k > val[new][i]) {
	      val[new][i] = val[last][j]-k;
	      dir[line][i] = x;
	    }
	    if (++j == lumres) j = 0;
	  }
	  k = offset[line] + i;
	  if (k >= lumres) k -= lumres;
	  /* silly edge detector: */
	  j = (+ lum[new][k - 1 + ((k-1 <       0) ? lumres : 0)]
	       + lum[new][k]
	       - lum[new][k + 1 - ((k+1 >= lumres) ? lumres : 0)]
	       - lum[new][k + 2 - ((k+2 >= lumres) ? lumres : 0)]);
	  if (j < 0) j = -j;
	  if (j > 10) j = 10 + (j-10)/4;
	  val[new][i] += j;
	  
	  /* don't allow the edge to get too near (< mindist) to a cutting
	   * point except in the first or last few lines */
	  if ((k < mindist || k >= lumres - mindist) && 
	      line >= 3 && line < lines - 3) {
	    val[new][i] -= 1000;
	    if (show_deadzone) rgb[line][k * resolution][0] = 255; 
	  }	  
	}
      }
    }
    if (verbose) 
      fprintf(stderr, "line %3d: offset = %3d (%.0f)\n",
	      linecount + ((halfframe == 2) ? 2 : 1), offset[line], maxsum);
    last ^= 1;  /* exchange buffers */
    new ^= 1;
    linecount += (halfframe ? 2 : 1);
  }
  
  if (corr_edge || mark_edge) {
    /* edge search, part 2 */
    maxsum = -1;
    for (i = 0; i < lumres; i++)   /* search starting point */
      if (val[last][i] > maxsum) x = i, maxsum = val[last][i];
    for (i = lines - 1; i >= 0; i--) {
      j = offset[i] + x;
      if (j > lumres) j -= lumres;
      if (j < 0) j += lumres;
      if (mark_edge) {
	rgb[i][j * resolution][0] = 0;
	rgb[i][j * resolution][1] = 255;
	rgb[i][j * resolution][2] = 0;
      }
      if (corr_edge) offset[i] = j;
      x += dir[i][x];
      if (x > lumres) x -= lumres;
      if (x < 0) x += lumres;
    }
  }
  
  /* write output image */
  if (verbose)
    fprintf(stderr, "writing data ...\n");
  for (line = 0; line < lines; line++) {
    if (markbreak) {
      hist[offset[line]]++;
      rgb[line][0][0] = rgb[line][0][1] = rgb[line][0][2] = 0;
      rgb[line][gout.xmax-1][0] = 255;
      rgb[line][gout.xmax-1][1] = 255;
      rgb[line][gout.xmax-1][2] = 255;
    }
    g_write_rgb(rgb[line] + lmargin, offset[line] * resolution, &gout);
    if (halfframe) 
      g_write_rgb(rgb[line] + lmargin, offset[line] * resolution, &gout);
  }
  fclose(gout.f);

  return 0;
}
