////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       Scientific pipeline I/O routine
//       Contributors:
//       Author: Andrea Zoli (IASF_Bologna)
//
// NOTICE
//       Any information contained in this software
//       is property of the AGILE TEAM and is strictly
//       private and confidential.
//       All rights reserved.
////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cmath>
#include <fitsio.h>
#include <pil.h>

using std::cout;
using std::endl;

void gaussian(double* kernel, int r) {
    int rr = 2*r+1;
    int ksz = rr*rr;
    double sigma = r/2.;

    double kt = 0;
    double a = 1./(sigma*sigma);
    double c = 1./(sigma*sigma);
    for (int y=-r; y<=r; y++) {
        for (int x=-r; x<=r; x++) {
            if ((x*x + y*y) <= r*r) {
                double v = exp(-.5*(a*x*x + c*y*y));
                kernel[(y+r)*rr+(x+r)] = v;
                kt += v;
            }
        }
    }

    // normalize kernel
    for (int aa=0; aa<ksz; aa++)
        kernel[aa] /= kt;
}

void convolve(double* kernel, double* src, double* dest,
              int width, int height, int r) {
    int rr = 2*r+1;

    double* dptr = dest;
    for (int jj=0; jj<height; jj++) {
        for (int ii=0; ii<width; ii++, dptr++) {
            for (int nn=jj-r, qq=0; nn<=jj+r; nn++, qq++) {
                if (nn>=0 && nn<height) {
                    register int nd = nn*width;
                    register int qd = qq*rr;
                    for (int mm=ii-r, pp=0; mm<=ii+r; mm++, pp++) {
                        if (mm>=0 && mm<width)
                            *dptr += src[nd+mm]*kernel[qd+pp];
                    }
                }
            }
        }
    }
}

void exitOnFitsErrors(int status) {
    if (status) {
        fits_report_error(stderr, status);
        exit(status);
    }
}

int main(int argc,char **argv) {
    int status = 0, numpar = 0;
    char mapFilename[FLEN_FILENAME];
    char kernFilename[FLEN_FILENAME];
    char outPrefix[FLEN_FILENAME];
    int smoothRadius = 0;

    status = PILInit(argc,argv);
    status = PILGetNumParameters(&numpar);
    status = PILGetString("mapfile", mapFilename);
    status = PILGetString("kernfile", kernFilename);
    status = PILGetInt("radius", &smoothRadius);
    status = PILGetString("outprefix", outPrefix);
    status = PILClose(status);

    cout << "#################################################################" << endl;
    cout << "########## AG_kernconv5.cpp v.1.0 - 20/01/2016          #########" << endl;
    cout << "#################################################################" << endl;
    cout << "#################################################################" << endl << endl;
    cout << "INPUT PARAMETERS:" << endl << endl;
    cout << "Map filename = " << mapFilename << endl;
    cout << "Kernel filename = " << kernFilename << endl;
    cout << "Output prefix = " << outPrefix << endl;
    cout << "Smoothing radius = " << smoothRadius << endl;

    cout << endl << "AG_kernconv...............................starting" << endl;

    // load the kernel matrix
    double simKernel[10000]; // max kernel 100x100
    std::ifstream is(kernFilename);
    int val = 0;
    char comma;
    double sumKern = 0.;
    for(val=0; is.good(); val++) {
        is >> simKernel[val] >> comma;
        sumKern += simKernel[val];
    }
    int kernSide = std::sqrt(val);
    int kernRadius = kernSide / 2;
    cout << "Kernel radius = " << kernRadius << endl;
    cout << "Kernel used = " << endl;
    for(int y=0; y<kernSide; y++) {
        for(int x=0; x<kernSide; x++) {
            simKernel[y*kernSide+x] /= sumKern;
            cout << simKernel[y*kernSide+x] << " ";
        }
        cout << endl;
    }
    cout << endl;

    // calculate source map with smooth
    double smoothKernel[10000]; // max kernel 100x100
    gaussian(smoothKernel, smoothRadius);
    int numaxis;
    long axis[2];
    fitsfile *infptr;
    fits_open_file(&infptr, mapFilename, READONLY, &status);
    exitOnFitsErrors(status);
    fits_get_img_dim(infptr, &numaxis, &status);
    exitOnFitsErrors(status);
    fits_get_img_size(infptr, 2, axis, &status);
    exitOnFitsErrors(status);
    int npixels = axis[0]*axis[1];
    double *img = new double[npixels];
    long pix1[3] = {1, 1, 1};
    fits_read_pix(infptr, TDOUBLE, pix1, npixels, NULL, img, NULL, &status);
    exitOnFitsErrors(status);
    double *imgSmooth = new double[npixels];
    convolve(simKernel, img, imgSmooth, axis[1], axis[0], kernRadius);
//    convolve(smoothKernel, img, imgSmooth, axis[1], axis[0], smoothRadius);
    bool saveMaps = true;
    if(saveMaps) {
        fitsfile *outfptr;
        std::string outfile(outPrefix);
        outfile += "_src_smooth.gz";
        fits_create_file(&outfptr, outfile.c_str(), &status);
        exitOnFitsErrors(status);
        fits_copy_header(infptr, outfptr, &status);
        exitOnFitsErrors(status);
        fits_write_pix(outfptr, TDOUBLE, pix1, npixels, imgSmooth, &status);
        fits_close_file(outfptr, &status);
        exitOnFitsErrors(status);
    }

    // get the peak of the smoothed source image within a fixed radius
    double peak = -1.;
    int centerX = axis[1]/2, centerY = axis[0]/2;
    int maxRadius = 50;
    int coordX = 0, coordY = 0;
    for (int y=0; y<axis[0]; ++y) {
        for (int x=0; x<axis[1]; ++x) {
            int index = y*axis[1] + x;
            if ((x - centerX)*(x - centerX) + (y - centerY)*(y - centerY) < maxRadius*maxRadius &&
                imgSmooth[index] > peak) {
                peak = imgSmooth[index];
                coordY = y;
                coordX = x;
            }
        }
    }

    // scale the sim kernel to the peak of the smoothed source image
    double simScaledKernel[10000];
    double scaleFactor = peak / simKernel[(kernSide/2)*kernSide+kernSide/2];
    double integralKernel = 0.0;
    double integralScaledKernel = 0.0;
    bool scale = true;
    if(scale) {
        for (int ky=0; ky<kernSide; ++ky) {
            for (int kx=0; kx<kernSide; ++kx) {
                int kindex = ky*kernSide + kx;
                simScaledKernel[kindex] = simKernel[kindex] * scaleFactor;
                integralKernel += simKernel[kindex];
                integralScaledKernel += simScaledKernel[kindex];
            }
        }
        cout << "Scaled sim kernel by a factor of " << scaleFactor << " is :" << endl;
        for (int ky=0; ky<kernSide; ++ky) {
            for (int kx=0; kx<kernSide; ++kx) {
                int kindex = ky*kernSide + kx;
                cout << simScaledKernel[kindex];
                cout << ", ";
            }
            cout << endl;
        }
        cout << "Integral of the sim kernel: " << integralKernel << endl;
        cout << "Integral of the scaled sim kernel: " << integralScaledKernel << endl;
    }
    else {
        for(int i=0; i<kernSide*kernSide; i++)
            simScaledKernel[i] = simKernel[i];
    }

    double *distImg = new double[npixels];
    for (int y=0; y<axis[0]; ++y) {
        for (int x=0; x<axis[1]; ++x) {
            int index = y*axis[1] + x;

            distImg[index] = 0.008;

            // skip border
            if ((x - centerX)*(x - centerX) + (y - centerY)*(y - centerY) >= maxRadius*maxRadius)
                continue;

            cout << "Source smooth img" << endl;
            cout << "y=" << y << " x=" << x << endl;
            for (int ky=0; ky<kernSide; ky++) {
                for (int kx=0; kx<kernSide; kx++) {
                    int dindex = (y-kernRadius+ky)*axis[1] + (x-kernRadius+kx);
                    cout << imgSmooth[dindex] << " ";
                }
                cout << endl;
            }

            distImg[index] = 0.0;

            // calculate the distance between the smoothed source image and
            // the simulated kernel scaled conveniently
            for (int ky=0; ky<kernSide; ky++) {
                for (int kx=0; kx<kernSide; kx++) {
                    int dindex = (y-kernRadius+ky)*axis[1] + (x-kernRadius+kx);
                    int kindex = ky*kernSide + kx;
                    distImg[index] += fabs(imgSmooth[dindex]-simScaledKernel[kindex]);
                    cout << simScaledKernel[kindex] << " " << imgSmooth[dindex] <<  " " << imgSmooth[dindex]-simScaledKernel[kindex] << endl;
/*                    distImg[index] += (imgSmooth[dindex]-simScaledKernel[kindex])*(imgSmooth[dindex]-simScaledKernel[kindex]);
                    cout << simScaledKernel[kindex] << " " << imgSmooth[dindex] <<  " " << (imgSmooth[dindex]-simScaledKernel[kindex])*(imgSmooth[dindex]-simScaledKernel[kindex]) << endl;*/
                }
            }
/*            distImg[index] = std::sqrt(distImg[index]); */
            cout << "dist total: " << distImg[index] << endl;
        }
    }
    if(saveMaps) {
        fitsfile *outfptr;
        std::string outfile(outPrefix);
        outfile += "_dist.gz";
        fits_create_file(&outfptr, outfile.c_str(), &status);
        exitOnFitsErrors(status);
        fits_copy_header(infptr, outfptr, &status);
        exitOnFitsErrors(status);
        fits_write_pix(outfptr, TDOUBLE, pix1, npixels, distImg, &status);
        fits_close_file(outfptr, &status);
        exitOnFitsErrors(status);
    }

    double minDist = 1.0;
    int minX = -1, minY = -1;
    for (int y=0; y<axis[0]; ++y) {
        for (int x=0; x<axis[1]; ++x) {
            int index = y*axis[1] + x;
            if(distImg[index] < minDist) {
                minDist = distImg[index];
                minX = x;
                minY = y;
            }
        }
    }
    std::ofstream os;
    std::string outfile(outPrefix);
    outfile += ".res";
    os.open(outfile.c_str(), std::ios_base::app);
    os << "minDist = " <<  minDist << " Y = " << minX << " X = " << minY << " peak = " << peak << " peakX = " << coordX << " peakY = " << coordY << endl;
    os.close();

    cout << "AG_kernconv............................... exiting" << endl;

    if(status) {
        cout << "AG_kernconv..................... exiting AG_kernconv ERROR:" << endl;
        fits_report_error(stdout, status);
    }
    else {
        cout << "\n\n\n#################################################################\n";
        cout <<       "#########  Task AG_kernconv........... exiting ##################\n";
        cout <<       "#################################################################\n\n\n";
    }

    return status;
}
