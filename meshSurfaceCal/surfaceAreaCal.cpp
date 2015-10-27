#include "objdata.h"

void fileHandle(char** dirpath);
int storeGeo(char* filepath, ObjGeo& objgeo);
void freeGeo(ObjGeo& objgeo);
float surfaceAreaCal(ObjGeo objgeo);
float dotProduct(Point PointA, Point PointB);
float mixedProduct(Point A, Point B, Point C);
Point crossProduct(Point PointA, Point PointB);

int main()
{
    char* filepath = new char[200];
    fileHandle(&filepath);
    ObjGeo objgeo;
    storeGeo(filepath, objgeo);
    float surfaceareavalue = 0;
    surfaceareavalue = surfaceAreaCal(objgeo);

    std::cout << "\n====================\nsurface area = " << surfaceareavalue  << "\n===================="<< std::endl;

    freeGeo(objgeo);
}

void fileHandle(char** filepath)
{
    std::cout << "Input obj filepath:" << std::endl;
    gets(*filepath);
}

int storeGeo(char* filepath, ObjGeo& objgeo)
{
    FILE* objfile;
    char* linetmp = new char[100]; *linetmp = '\0';
    int Vcount = 0, Fcount = 0, VNcount = 0;
    if (0 != (objfile = fopen(filepath, "r")))
    {
        while (0 != fgets(linetmp, 100, objfile)) //count Vnum and Fnum of every parts, no more than 200 organs
        {
            if (1 == strlen(linetmp))
                continue;
            *(linetmp + strlen(linetmp) - 1) = '\0';
            if ('v' == *linetmp && ' ' == *(linetmp + 1))
            {
                Vcount++;
            }
            else if ('v' == *linetmp && 'n' == *(linetmp + 1))
            {
                VNcount++;
            }
            else if ('f' == *linetmp && ' ' == *(linetmp + 1))
            {
                Fcount++;
            }
        }

        //malloc enough space
        objgeo.facecount = Fcount;
        objgeo.pointcount = Vcount;
        objgeo.vncount = VNcount;
        objgeo.point = new Point[Vcount + 1];
        objgeo.pointnormal = new Point[VNcount + 1];
        objgeo.face = new Face[Fcount + 1];

        //back to the begning of file
        fseek(objfile, 0, SEEK_SET);
        int pointtmp = 0; int pointnormaltmp = 0; int facetmp = 0;
        while (0 != fgets(linetmp, 100, objfile)) //count Vnum and Fnum of every parts, no more than 200 organs
        {
            if (1 == strlen(linetmp))
                continue;
            *(linetmp + strlen(linetmp) - 1) = '\0';
            if ('v' == *linetmp && ' ' == *(linetmp + 1))
            {
                float coorX, coorY, coorZ;
                sscanf(linetmp, "%*s %f%f%f", &coorX, &coorY, &coorZ);
                (objgeo.point + pointtmp + 1)->setcoorX(coorX);
                (objgeo.point + pointtmp + 1)->setcoorY(coorY);
                (objgeo.point + pointtmp + 1)->setcoorZ(coorZ);
                pointtmp++;
            }
            else if ('v' == *linetmp && 'n' == *(linetmp + 1))
            {
                float coorX, coorY, coorZ;
                sscanf(linetmp, "%*s %f%f%f", &coorX, &coorY, &coorZ);
                (objgeo.pointnormal + pointnormaltmp + 1)->setcoorX(coorX);
                (objgeo.pointnormal + pointnormaltmp + 1)->setcoorY(coorY);
                (objgeo.pointnormal + pointnormaltmp + 1)->setcoorZ(coorZ);
                pointnormaltmp++;
            }
            else if ('f' == *linetmp && ' ' == *(linetmp + 1))
            {
                int facetmp1, facetmp2, facetmp3;
                sscanf(linetmp, "%*s %d//%*d %d//%*d %d//%*d", &facetmp1, &facetmp2, &facetmp3);
                (objgeo.face + facetmp)->setpointAindex(facetmp1);
                (objgeo.face + facetmp)->setpointBindex(facetmp2);
                (objgeo.face + facetmp)->setpointCindex(facetmp3);
                facetmp++;
            }
        }

        fclose(objfile);
        return 0;
    }
    else
    {
        return 1;
    }
}

void freeGeo(ObjGeo& objgeo)
{
    delete[]objgeo.face;
    delete[]objgeo.point;
    delete[]objgeo.pointnormal;
}

// Calculate total surface area value of triangles that construct our mesh
float surfaceAreaCal(ObjGeo objgeo)
{
    float surfaceareavalue = 0;
    for (int i = 0; i < objgeo.facecount; i++)
    {
        int Aindex = (objgeo.face + i)->getpointAindex();
        int Bindex = (objgeo.face + i)->getpointBindex();
        int Cindex = (objgeo.face + i)->getpointCindex();
        // We should use cross product operation instead of dot product operation
        Point BA, BC, TMP;
        BA.setcoorX((objgeo.point + Aindex)->getcoorX() - (objgeo.point + Bindex)->getcoorX());
        BA.setcoorY((objgeo.point + Aindex)->getcoorY() - (objgeo.point + Bindex)->getcoorY());
        BA.setcoorZ((objgeo.point + Aindex)->getcoorZ() - (objgeo.point + Bindex)->getcoorZ());

        BC.setcoorX((objgeo.point + Cindex)->getcoorX() - (objgeo.point + Bindex)->getcoorX());
        BC.setcoorY((objgeo.point + Cindex)->getcoorY() - (objgeo.point + Bindex)->getcoorY());
        BC.setcoorZ((objgeo.point + Cindex)->getcoorZ() - (objgeo.point + Bindex)->getcoorZ());

        TMP = crossProduct(BA, BC);

        surfaceareavalue += 0.5 * sqrt(dotProduct(TMP, TMP));

    }

    return surfaceareavalue;
}

// Calculate dot Product
// Applicable to 3 dims dot product
float dotProduct(Point PointA, Point PointB)
{
    float dotproductans = 0;
    dotproductans += PointA.getcoorX() * PointB.getcoorX();
    dotproductans += PointA.getcoorY() * PointB.getcoorY();
    dotproductans += PointA.getcoorZ() * PointB.getcoorZ();

    return dotproductans;
}

//Calculate (A x B)*C
float mixedProduct(Point A, Point B, Point C)
{
    Point D; //A x B = D
    D.setcoorX(A.getcoorY() * B.getcoorZ() - A.getcoorZ() * B.getcoorY());
    D.setcoorY(A.getcoorZ() * B.getcoorX() - A.getcoorX() * B.getcoorZ());
    D.setcoorZ(A.getcoorX() * B.getcoorY() - A.getcoorY() * B.getcoorX());

    return (D.getcoorX() * C.getcoorX() + D.getcoorY() * C.getcoorY() + D.getcoorZ() * C.getcoorZ());
}

//Calculat PointA x PointB
Point crossProduct(Point PointA, Point PointB)
{
    Point pointans;
    pointans.setcoorX(PointA.getcoorY() * PointB.getcoorZ() - PointB.getcoorY() * PointA.getcoorZ());
    pointans.setcoorY(PointA.getcoorZ() * PointB.getcoorX() - PointB.getcoorZ() * PointA.getcoorX());
    pointans.setcoorZ(PointA.getcoorX() * PointB.getcoorY() - PointB.getcoorX() * PointA.getcoorY());

    return pointans;
}