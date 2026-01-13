set -v
P=$1
m=$2
sed -i 's/#include <filesystem>//g' $P.cpp
sed -i 's/std::filesystem::path(dataPath) \/ \"/\"learned-systems-rmi\/rmi_data\//g' $P.cpp
sed -i '/void cleanup/i bool save(char const *filename);' $P.h

sed -i "/void cleanup/i \\
bool save(char const *filename) {\\
    int64_t L1_SIZE = ${m};\\
    std::ofstream outstream(filename, std::ofstream::binary);\\
    outstream.seekp(0);\\
    outstream.write((char*)&L0_PARAMETER0, sizeof(double));\\
    outstream.write((char*)&L0_PARAMETER1, sizeof(double));\\
    outstream.write((char*)&L0_PARAMETER2, sizeof(double));\\
    outstream.write((char*)&L0_PARAMETER3, sizeof(double));\\
    outstream.write((char*)&L0_PARAMETER4, sizeof(double));\\
    outstream.write((char*)&L0_PARAMETER5, sizeof(double));\\
    outstream.write((char*)&L0_PARAMETER6, sizeof(double));\\
    outstream.write((char*)&L0_PARAMETER7, sizeof(double));\\
    outstream.write((char*)&L0_PARAMETER8, sizeof(double));\\
    outstream.write((char*)&L0_PARAMETER9, sizeof(double));\\
    outstream.write((char*)&L0_PARAMETER10, sizeof(double));\\
    outstream.write((char*)&L0_PARAMETER11, sizeof(double));\\
    outstream.write((char*)&L0_PARAMETER12, sizeof(double));\\
    outstream.write((char*)&L0_PARAMETER13, sizeof(double));\\
    outstream.write((char*)&L0_PARAMETER14, sizeof(double));\\
    outstream.write((char*)&L0_PARAMETER15, sizeof(double));\\
    outstream.write((char*)&L0_PARAMETER16, sizeof(double));\\
    outstream.write((char*)&L0_PARAMETER17, sizeof(double));\\
    outstream.write((char*)&L0_PARAMETER18, sizeof(double));\\
    outstream.write((char*)&L0_PARAMETER19, sizeof(double));\\
    outstream.write((char*)&L0_PARAMETER20, sizeof(double));\\
    outstream.write((char*)&L0_PARAMETER21, sizeof(double));\\
    outstream.write((char*)&L0_PARAMETER22, sizeof(double));\\
    outstream.write((char*)&L1_SIZE, sizeof(int64_t));\\
    outstream.write((char*)L1_PARAMETERS, L1_SIZE * 3 * sizeof(double));\\
    outstream.close();\\
    return true;\\
}\\
" $P.cpp