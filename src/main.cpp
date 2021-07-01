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
#include <chrono>

typedef uint32_t Color;
typedef uint32_t TileIndex;
TileIndex TileMask{0xFFFFFFFF};

uint64_t qwordPerTiles{};

struct Parameters
{
    uint16_t tileWidth, tileHeight;
    std::string filename;
    uint16_t mapWidth, mapHeight;
};

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
size_t tilesSize;
struct Tile
{
    std::vector<Color> bitmap;
    //std::vector<bool> compatibility[4]; // +x, +y, -x, -y
    uint64_t* compatibilityFast;

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

bool ParseParameters(int argc, char** argv, Parameters& parameters)
{
    parameters.tileWidth = 16;
    parameters.tileHeight = 16;

    parameters.mapWidth = 26;
    parameters.mapHeight = 26;

#ifdef _MSC_VER
    parameters.filename = "../maps/Zelda3LightOverworldBG_masked.png";
    //parameters.filename = "../maps/Zelda3LightOverworldBG.png";
#else
    char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd)) != NULL)
    {
        printf("Current working dir: %s\n", cwd);
    }
    parameters.filename = std::string(cwd) + "/../../maps/Zelda3LightOverworldBG_masked.png";
#endif
    return true;
}

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
    int res = -1;
    for (int i = 0; i < tileCount; i++)
    {
        //if (mCoef[idx + i])
        if (GetCoef(coord, i))
        {
            assert(res == -1);
            res = i;
        }
    }
    return res;
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
    //int idx = (coord.y * mWidth + coord.x) * tileCount;
    for (uint64_t i = 0; i < tileCount; i++)
    {
        //if (mCoef[idx + i])
        if (GetCoef(coord, i))
        {
            mTotalSum--;
        }
        //mCoef[idx + i] = 0;
        SetCoef0(coord, i);
    }
    //mCoef[idx + tileIndex] = 1;
    SetCoef1(coord, tileIndex);
    mTotalSum++;
    mSumCoef[coord.y * mWidth + coord.x] = 1;
}
void Collapse(Coord coord)
{
    int *potentials = new int [tileCount];
    int potentialIndex = 0;
    //int idx = (coord.y * mWidth + coord.x) * tileCount;
    //auto & v = mCoef[coord.y * mWidth + coord.x];
    int cnt = 0;
    for (uint64_t i = 0; i < tileCount; i++)
    {
        //if (mCoef[idx + i])
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
    //int idx = (coord.y * mWidth + coord.x) * tileCount;
    //auto & v = mCoef[coord.y * mWidth + coord.x];
    //assert(mCoef[idx + tileIndex]);
    //mCoef[idx + tileIndex] = 0;
    assert(GetCoef(coord, tileIndex));
    SetCoef0(coord, tileIndex);
    mSumCoef[coord.y * mWidth + coord.x] --;
    mTotalSum--;
}

bool TileCompatibleFast(const Tile& tile1, TileIndex tileIndex2, int direction)
{
    //int key1 = mTiles[tileIndex1].mKeys[GetAngle(dir)];
    //int key2 = mTiles[tileIndex2].mKeys[GetHook(dir)];

    //return (key1 & key2) != 0;
    return tile1.GetCompatibilityFast(direction, tileIndex2);
}


int GetPossibleTiles(Coord coord, int* possibleTiles)
{
    int res = 0;
    //int idx = (coord.y * mWidth + coord.x) * tileCount;
    //auto & v = mCoef[coord.y * mWidth + coord.x];
    for (int i = 0; i < tileCount; i++)
    {
        //if (mCoef[idx + i])
        if (GetCoef(coord, i))
        {
            possibleTiles[res++] = i;
        }
    }
    assert(res);
    return res;
}

struct BumpAllocator
{
    BumpAllocator()
    {
        m_ptrSource = m_ptr = (uint8_t*)malloc(PageAllocation);
        m_totalAllocated = 0;
    }
    
    template<typename T> T* Allocate(size_t itemCount)
    {
        return (T*)DoAllocate(itemCount * sizeof(T));
    }
    
    void* DoAllocate(size_t sizeInBytes)
    {
        m_totalAllocated += sizeInBytes;
        assert(m_totalAllocated < PageAllocation);
        void* res = m_ptr;
        m_ptr += sizeInBytes;
        return res;
    }
    
    void Reset()
    {
        m_ptr = m_ptrSource;
        m_totalAllocated = 0;
    }
    uint8_t* m_ptr;
    uint8_t* m_ptrSource;
    size_t m_totalAllocated;
    static const size_t PageAllocation = 100 * 1024 * 1024;
};

void Propagate(Coord coord)
{
    static std::vector<Coord> coords;
    coords.clear();
    coords.push_back(coord);
    static BumpAllocator bump;
    int infoTick = 0;
    while (coords.size())
    {
        infoTick ++;
        if (infoTick >= 1000)
        {
            printf("Remaining : %d Entropy : %d     \r", int(coords.size()), int(mTotalSum));
            infoTick = 0;
        }
        bump.Reset();
        Coord currentCoord = coords.back();
        coords.pop_back();

        int* curPossibleTiles = bump.Allocate<int>(tileCount);
        int curPossibleTileCount = GetPossibleTiles(currentCoord, curPossibleTiles);
        uint64_t* currentCoefs = GetCoefs(currentCoord);
        Coord validDirs[4];
        int directions[4];
        int validDirCount = GetValidDirs(currentCoord, validDirs, directions, mWidth, mHeight);
        for (int d = 0; d < validDirCount; d++)
        {
            static const int inverser[4] = { 2,3,0,1 };
            const auto inversedDir = inverser[directions[d]];
            Coord dir = validDirs[d];
            Coord otherCoord = { currentCoord.x + dir.x, currentCoord.y + dir.y };
            int* otherPossibleTiles = bump.Allocate<int>(tileCount);
            int otherPossibleTileCount = GetPossibleTiles(otherCoord, otherPossibleTiles);
            for (int otherTileIndex = 0; otherTileIndex < otherPossibleTileCount; otherTileIndex++)
            {
                int otherTile = otherPossibleTiles[otherTileIndex];
                const Tile& tile1 = tiles[otherTile];
                const uint64_t* otherCoefs = tile1.GetCompatibilityCoefs(inversedDir);
                bool tileCompatible = false;

                for (size_t batch = 0; batch < qwordPerTiles; batch++)
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

void Blit(const Color* tileInput, Color* output, uint32_t offset, uint32_t pitch, uint16_t tileWidth, uint16_t tileHeight)
{
    for(size_t y = 0; y < tileHeight; y++)
    {
        for(size_t x = 0; x < tileWidth; x++)
        {
            Color* destination = output + pitch * y + offset + x;
            *destination = *tileInput++;
        }
    }
}

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
    tilesSize = tiles.size();
    //tilesPtr = tiles.data();
    
    // resize compatibility list
    qwordPerTiles = (tiles.size() + 63ULL) >> 6ULL;
    const uint64_t qwordPerTilesAllDirections = qwordPerTiles * 4ULL;
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
    //mCoef.resize(mWidth * mHeight * tileCount, true);
    mCoefxx = (uint64_t*)malloc(coefTotalBytes);


    for (int mapIndex = 0; mapIndex < 10; mapIndex++)
    {

        for (size_t i = 0; i < totalQWords; i++)
        {
            mCoefxx[i] = 0xFFFFFFFFFFFFFFFFULL;
        }
        mSumCoef.clear();
        mSumCoef.resize(mWidth * mHeight, tileCount);
        mTotalSum = mWidth * mHeight * tileCount;

        srand(7888 + mapIndex);

        try
        {
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
            std::vector<Color> image(mWidth * parameters.tileWidth * mHeight * parameters.tileHeight);

            for (size_t ty = 0; ty < mHeight; ty++)
            {
                for (size_t tx = 0; tx < mWidth; tx++)
                {
                    Coord coord{ int32_t(tx), int32_t(ty) };
                    const auto tileIndex = GetTileAtIndex(coord);
                    const auto& tile = tiles[tileIndex];
                    const uint32_t offset = uint32_t(ty * parameters.tileHeight * mWidth * parameters.tileWidth + tx * parameters.tileWidth);
                    Blit(tile.bitmap.data(), image.data(), offset, mWidth * parameters.tileWidth, parameters.tileWidth, parameters.tileHeight);
                }
            }

            // save image
            const auto imageWidth = mWidth * parameters.tileWidth;
            char mapFilename[512];
            sprintf(mapFilename, "MapGen_%d.png", mapIndex);
#ifdef _MSC_VER
            int res = stbi_write_png(mapFilename, imageWidth, mHeight * parameters.tileHeight, 4,
                image.data(), imageWidth * 4);
#else
            int res = stbi_write_png("/Users/cedricguillemet/dev/WFCTiles/res.png", imageWidth, mHeight * parameters.tileHeight, 4,
                image.data(), imageWidth * 4);
#endif
        }
        catch(...)
        {
            printf("Failing at generating %d\n", mapIndex);
        }
    }
    return 0;
}
