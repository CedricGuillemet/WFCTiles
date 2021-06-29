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

typedef uint32_t Color;
typedef uint32_t TileIndex;
TileIndex TileMask{0xFFFFFFFF};

struct Parameters
{
    uint16_t tileWidth, tileHeight;
    std::string filename;
};

struct Coord
{
    int32_t x, y;
};

int GetValidDirs(Coord coord, Coord *dest, int* directions, uint32_t tileXCount, uint32_t tileYCount)
{
    int coordCount = 0;
    if (coord.x < (tileXCount-1))
    {
        directions[coordCount] = 0;
        dest[coordCount++] = {1, 0 };
    }
    if (coord.y < (tileYCount - 1))
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
    std::vector<bool> compatibility[4]; // +x, +y, -x, -y
    
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
#ifdef _MSC_VER
    parameters.filename = "../maps/Zelda3LightOverworldBG_masked.png";
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
    const size_t tileXCount = x / parameters.tileWidth;
    const size_t tileYCount = y / parameters.tileHeight;
    tileMap.resize(tileXCount * tileYCount);

    for (size_t ty = 0; ty < tileYCount; ty++)
    {
        for (size_t tx = 0; tx < tileYCount; tx++)
        {
            TileIndex tileIndex = GetTileIndex((Color*)data, (ty * parameters.tileHeight) * x + tx * parameters.tileWidth, x, parameters);
            tileMap[ty * tileXCount + tx] = tileIndex;
        }
    }
    printf("%d x %d source tiles with %d unique tiles\n", int(tileXCount), int(tileYCount), int(tiles.size()));
    
    /*
    stbi_write_png("tile 0.png", 16, 16, 4, tiles[0].bitmap.data(), 16 * 4);
    stbi_write_png("tile 1.png", 16, 16, 4, tiles[1].bitmap.data(), 16 * 4);
    stbi_write_png("tile 1516.png", 16, 16, 4, tiles[1516].bitmap.data(), 16 * 4);
    stbi_write_png("tile 1517.png", 16, 16, 4, tiles[1517].bitmap.data(), 16 * 4);
     */
    // resize compatibility list
    for (auto& tile : tiles)
    {
        for (auto& compatibility : tile.compatibility)
        {
            compatibility.resize(tiles.size(), false);
        }
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
                tiles[tileIndex].compatibility[validDirectionIndex[i]][neighboorTileIndex] = true;
            }
        }
    }
    // check compatibility
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

    stbi_image_free(data);
    return 0;
}
