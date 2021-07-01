#include <stdio.h>
#ifndef _MSC_VER
#include <unistd.h>
#endif
#include <string>
#include <vector>
#include <map>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "Parameters.h"
#include <chrono>

typedef uint32_t Color;
typedef uint32_t TileIndex;
TileIndex TileMask{0xFFFFFFFF};

struct Coord
{
    union {
        struct {
            int32_t x, y;
        };
        uint64_t val;
    };

    const bool operator ==(const Coord& other) const
    {
        return val == other.val;
    }
};

int GetValidDirs(Coord coord, Coord *dest, int* directions, uint32_t tileXCount, uint32_t tileYCount)
{
    int coordCount = 0;
    if (coord.x < int32_t(tileXCount-1))
    {
        directions[coordCount] = 0;
        dest[coordCount++] = {1, 0 };
    }
    if (coord.y < int32_t(tileYCount - 1))
    {
        directions[coordCount] = 1;
        dest[coordCount++] = { 0, 1 };
    }
    if (coord.x > 0)
    {
        directions[coordCount] = 2;
        dest[coordCount++] = { -1, 0 };
    }
    if (coord.y > 0)
    {
        directions[coordCount] = 3;
        dest[coordCount++] = {0, -1};
    }
    return coordCount;
}

struct Tile
{
    std::vector<Color> bitmap;
    uint64_t* compatibilityFast;

    static inline uint64_t qwordPerTiles{};

    const uint64_t* GetCompatibilityCoefs(size_t direction) const
    {
        return &compatibilityFast[direction * qwordPerTiles];
    }

    void SetCompatibilityFast(size_t direction, size_t index)
    {
        uint64_t bitIndex = direction * qwordPerTiles * 64ULL + index;
        uint64_t byteIndex = bitIndex >> 6ULL;
        uint64_t bitInByte = bitIndex & 63ULL;
        compatibilityFast[byteIndex] |= 1ULL << bitInByte;
    }

    bool GetCompatibilityFast(size_t direction, size_t index) const
    {
        uint64_t bitIndex = direction * qwordPerTiles * 64ULL + index;
        uint64_t byteIndex = bitIndex >> 6ULL;
        uint64_t bitInByte = bitIndex & 63ULL;
        return compatibilityFast[byteIndex] & (1ULL << bitInByte);
    }

    /*
    bool HasCompatible(int direction) const
    {
        const auto& comp = compatibility[direction];
        for(size_t i = 0; i < comp.size(); i++)
        {
            if (comp[i])
            {
                return true;
            }
        }
        return false;
    }
    */
    uint64_t ComputeChecksum() const
    {
        uint64_t res = 0;
        for(const auto color : bitmap)
        {
            res += color;
        }
        return res;
    }
    
    bool IsMask() const
    {
        for (size_t i = 0; i < bitmap.size(); i++)
        {
            if (bitmap[i] != 0xFFFF00FF)
            {
                return false;
            }
        }
        return true;
    }
    bool operator == (const Tile& other) const
    {
        assert(bitmap.size() == other.bitmap.size());
        for(size_t i = 0; i < bitmap.size(); i++)
        {
            if (bitmap[i] != other.bitmap[i])
            {
                return false;
            }
        }
        return true;
    }
    
    void BuildFromMap(Color* colors, uint32_t offset, uint32_t pitch, const Parameters& parameters)
    {
        const auto w = parameters.tileWidth;
        const auto h = parameters.tileHeight;
        bitmap.resize(h * w);
        auto* p = bitmap.data();
        for (size_t y = 0; y < h; y++)
        {
            for (size_t x = 0; x < w; x++)
            {
                *p++ = colors[offset + y * pitch +x];
            }
        }
    }
};

std::vector<Tile> tiles;
std::map<uint64_t, std::vector<TileIndex> > uniqueTiles;

TileIndex GetTileIndex(Color* colors, uint32_t offset, uint32_t pitch, const Parameters& parameters)
{
    Tile tile;
    tile.BuildFromMap(colors, offset, pitch, parameters);
    if (tile.IsMask())
    {
        return TileMask;
    }
    auto checkSum = tile.ComputeChecksum();
    
    auto iter = uniqueTiles.find(checkSum);
    if (iter != uniqueTiles.end())
    {
        auto& existingTileIndices = iter->second;
        for (const auto tileIndex : existingTileIndices)
        {
            if (tiles[tileIndex] == tile)
            {
                // tile found
                return tileIndex;
            }
        }
    }
    // not found, add new
    TileIndex res = static_cast<TileIndex>(tiles.size());
    tiles.push_back(tile);
    uniqueTiles[checkSum].push_back(res);
    return res;
}

std::vector<TileIndex> tileMap;

int mWidth, mHeight;
uint64_t* mCoefxx{};
size_t coefQwordsPerTile;
std::vector<unsigned short> mSumCoef;
unsigned int mTotalSum;
int tileCount;

bool GetCoef(Coord coord, uint64_t index)
{
    const uint64_t bitIndex = (coord.y * mWidth + coord.x) * coefQwordsPerTile * 64 + index;
    const uint64_t qwordIndex = bitIndex >> 6ULL;
    const uint64_t bitInQword = bitIndex & 63ULL;
    return mCoefxx[qwordIndex] & (1ULL << bitInQword);
}

uint64_t* GetCoefs(Coord coord)
{
    return &mCoefxx[(coord.y * mWidth + coord.x) * coefQwordsPerTile];
}

void SetCoef0(Coord coord, uint64_t index)
{
    const uint64_t bitIndex = (coord.y * mWidth + coord.x) * coefQwordsPerTile * 64 + index;
    const uint64_t qwordIndex = bitIndex >> 6ULL;
    const uint64_t bitInQword = bitIndex & 63ULL;
    mCoefxx[qwordIndex] &= ~(1ULL << bitInQword);
}

void SetCoef1(Coord coord, uint64_t index)
{
    const uint64_t bitIndex = (coord.y * mWidth + coord.x) * coefQwordsPerTile * 64 + index;
    const uint64_t qwordIndex = bitIndex >> 6ULL;
    const uint64_t bitInQword = bitIndex & 63ULL;
    mCoefxx[qwordIndex] |= (1ULL << bitInQword);
}

int GetTileAtIndex(Coord coord)
{
    for (int i = 0; i < tileCount; i++)
    {
        if (GetCoef(coord, i))
        {
            return i;
        }
    }
    return -1;
}

Coord GetMinEntropy()
{
    int minEntropy = INT_MAX;
    static std::vector<Coord> minEntropyCoords;
    minEntropyCoords.clear();
    for (int y = 0; y < mHeight; y++)
    {
        for (int x = 0; x < mWidth; x++)
        {
            int coef = mSumCoef[y * mWidth + x];
            if (coef == 1)
                continue;
            if (coef < minEntropy)
            {
                minEntropy = coef;
                minEntropyCoords.clear();
            }
            if (coef == minEntropy)
                minEntropyCoords.push_back({ x, y });
        }
    }
    assert(!minEntropyCoords.empty());
    return minEntropyCoords[rand() % minEntropyCoords.size()];
}
void Collapse(Coord coord, int tileIndex)
{
    for (uint64_t i = 0; i < tileCount; i++)
    {
        if (GetCoef(coord, i))
        {
            mTotalSum--;
        }
        SetCoef0(coord, i);
    }
    SetCoef1(coord, tileIndex);
    mTotalSum++;
    mSumCoef[coord.y * mWidth + coord.x] = 1;
}
void Collapse(Coord coord)
{
    int *potentials = new int [tileCount];
    int potentialIndex = 0;

    int cnt = 0;
    for (uint64_t i = 0; i < tileCount; i++)
    {
        if (GetCoef(coord, i))
        {
            potentials[potentialIndex++] = int(i);
        }
    }
    assert(potentialIndex);
    static int rd = 0;
    int selected = potentials[rand() % potentialIndex];
    delete [] potentials;
    Collapse(coord, selected);
}
bool IsFullyCollapsed()
{
    return mTotalSum == (mWidth * mHeight);
}
void Constrain(Coord coord, int tileIndex)
{
    assert(GetCoef(coord, tileIndex));
    SetCoef0(coord, tileIndex);
    mSumCoef[coord.y * mWidth + coord.x] --;
    mTotalSum--;
}

int GetPossibleTiles(Coord coord, int* possibleTiles)
{
    int res = 0;
    for (int i = 0; i < tileCount; i++)
    {
        if (GetCoef(coord, i))
        {
            possibleTiles[res++] = i;
        }
    }
    assert(res);
    return res;
}

void Propagate(Coord coord)
{
    static std::vector<Coord> coords;
    coords.clear();
    coords.push_back(coord);
    int infoTick = 0;
    static int* curPossibleTiles{ new int[tileCount] };
    while (coords.size())
    {
        infoTick ++;
        if (infoTick >= 1000)
        {
            printf("Remaining : %d Entropy : %d     \r", int(coords.size()), int(mTotalSum));
            infoTick = 0;
        }
        Coord currentCoord = coords.back();
        coords.pop_back();
        
        int curPossibleTileCount = GetPossibleTiles(currentCoord, curPossibleTiles);
        uint64_t* currentCoefs = GetCoefs(currentCoord);
        Coord validDirs[4];
        int directions[4];
        int validDirCount = GetValidDirs(currentCoord, validDirs, directions, mWidth, mHeight);
        for (int d = 0; d < validDirCount; d++)
        {
            static const int inverser[4] = { 2,3,0,1 };
            const auto inversedDir = inverser[directions[d]];
            const Coord dir = validDirs[d];
            const Coord otherCoord = { currentCoord.x + dir.x, currentCoord.y + dir.y };
            
            const uint64_t* otherMapCoefs = GetCoefs(otherCoord);
            
            for (TileIndex otherTile = 0; otherTile < TileIndex(tileCount); otherTile++)
            {
                uint64_t otherQword = otherTile >> 6ULL;
                uint64_t otherBit = 1ULL << (otherTile & 63ULL);
                if ((otherMapCoefs[otherQword] & otherBit) == 0)
                {
                    continue;
                }

                const Tile& tile1 = tiles[otherTile];

                const uint64_t* otherCoefs = tile1.GetCompatibilityCoefs(inversedDir);
                bool tileCompatible = false;

                for (size_t batch = 0; batch < Tile::qwordPerTiles; batch++)
                {
                    if (otherCoefs[batch] & currentCoefs[batch])
                    {
                        tileCompatible = true;
                        break;
                    }
                }
                if (!tileCompatible)
                {
                    Constrain(otherCoord, otherTile);
                    {
                        auto iter = std::find(coords.begin(), coords.end(), otherCoord);
                        if (iter == coords.end())
                        {
                            coords.push_back(otherCoord);
                        }
                    }
                }
            }
        }
    }
    //printf("\nMax bump %d\n", int(bump.m_totalAllocated));
}


#include "MapImage.h"

int main(int argc, char** argv)
{
    printf("WFCTiles\n");
    
    Parameters parameters;
    if (!ParseParameters(argc, argv, parameters))
    {
        return 1;
    }

    int x,y,n;
    unsigned char *data = stbi_load(parameters.filename.c_str(), &x, &y, &n, 4);
    if (!data)
    {
        return 1;
    }

    // build tile map
    const uint32_t tileXCount = x / parameters.tileWidth;
    const uint32_t tileYCount = y / parameters.tileHeight;
    tileMap.resize(tileXCount * tileYCount);

    for (uint32_t ty = 0; ty < tileYCount; ty++)
    {
        for (uint32_t tx = 0; tx < tileYCount; tx++)
        {
            uint32_t offset = (ty * parameters.tileHeight) * x + tx * parameters.tileWidth;
            TileIndex tileIndex = GetTileIndex((Color*)data, offset, x, parameters);
            tileMap[ty * tileXCount + tx] = tileIndex;
        }
    }
    printf("%d x %d source tiles with %d unique tiles\n", int(tileXCount), int(tileYCount), int(tiles.size()));
    stbi_image_free(data);
    
    // resize compatibility list
    Tile::qwordPerTiles = (tiles.size() + 63ULL) >> 6ULL;
    const uint64_t qwordPerTilesAllDirections = Tile::qwordPerTiles * 4ULL;
    const uint64_t qwordPerTilesAllDirectionsInBytes = qwordPerTilesAllDirections * sizeof(uint64_t);
    for (auto& tile : tiles)
    {
        tile.compatibilityFast = (uint64_t*)malloc(qwordPerTilesAllDirectionsInBytes);
        memset(tile.compatibilityFast, 0, qwordPerTilesAllDirectionsInBytes);
    }
    
    // build compatibility
    for (size_t ty = 0; ty < tileYCount; ty++)
    {
        for (size_t tx = 0; tx < tileXCount; tx++)
        {
            TileIndex tileIndex = tileMap[ty * tileXCount + tx];
            if (tileIndex == TileMask)
            {
                continue;
            }
            Coord coord{int(tx), int(ty)};
            Coord validDirections[4];
            int validDirectionIndex[4];
            int directionCount = GetValidDirs(coord, validDirections, validDirectionIndex, tileXCount, tileYCount);
            for (int i = 0; i < directionCount; i++)
            {
                Coord neighbourCoord = coord;
                neighbourCoord.x += validDirections[i].x;
                neighbourCoord.y += validDirections[i].y;
                TileIndex neighboorTileIndex = tileMap[neighbourCoord.y * tileXCount + neighbourCoord.x];
                if (neighboorTileIndex == TileMask)
                {
                    continue;
                }
                //tiles[tileIndex].compatibility[validDirectionIndex[i]][neighboorTileIndex] = true;
                tiles[tileIndex].SetCompatibilityFast(validDirectionIndex[i], neighboorTileIndex);
                
            }
        }
    }
    
    // check compatibility
    /*
    printf("Checking compatibility ...\n");
    for (size_t i = 0; i<tiles.size(); i++)
    {
        for (int direction = 0; direction < 4; direction++)
        {
            bool hasCompatible = tiles[i].HasCompatible(direction);
            if (!hasCompatible)
            {
                printf("Tile %d has no compatibility for direction %d\n", int(i), direction);
            }
        }
    }
    printf("Done\n");
    */
    // do it
    mWidth = parameters.mapWidth;
    mHeight = parameters.mapHeight;
    tileCount = int(tiles.size());

    const size_t coefBitsPerTile = (tileCount + 64) & (~63);
    coefQwordsPerTile = coefBitsPerTile >> 6;
    const size_t totalQWords = mWidth * mHeight * coefQwordsPerTile;
    const size_t coefTotalBytes = totalQWords * sizeof(uint64_t);
    mCoefxx = (uint64_t*)malloc(coefTotalBytes);

    for (int mapIndex = 0; mapIndex < 1; mapIndex++)
    {
        for (size_t i = 0; i < totalQWords; i++)
        {
            mCoefxx[i] = 0xFFFFFFFFFFFFFFFFULL;
        }
        mSumCoef.clear();
        mSumCoef.resize(mWidth * mHeight, tileCount);
        mTotalSum = mWidth * mHeight * tileCount;

        srand(parameters.seed + mapIndex);

            auto startTime = std::chrono::high_resolution_clock::now();
            while (!IsFullyCollapsed())
            {
                Coord coord = GetMinEntropy();
                Collapse(coord);
                Propagate(coord);
            }
            auto endTime = std::chrono::high_resolution_clock::now();
            auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);

            printf("Time for collapsing all: %2.4f seconds\n", float(time_span.count()));

            // generate image
            SaveMapImage(parameters, mapIndex);
    }
    return 0;
}
