#include <string>
#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <map>
#include <vector>
#include <array>
#include <queue>
#include <sstream>
#include <iomanip>
#include <math.h>
#include "needlemanWunsch.cpp"
#include "initData.cpp"
using namespace std;

void swapDouble(double* a, double* b)
{
    double t = *a;
    *a = *b;
    *b = t;
}

void swapString(string* a, string* b)
{
    string t = *a;
    *a = *b;
    *b = t;
}

int partition (vector<string>& names, vector<double>& data, int low, int high)
{
    double pivot = data[high];
    int i = (low - 1);
 
    for (int j = low; j <= high- 1; j++)
    {
        if (data[j] <= pivot)
        {
            i++;   
            swapDouble(&data[i], &data[j]);
            swapString(&names[i], &names[j]);
        }
    }
    swapDouble(&data[i + 1], &data[high]);
    swapString(&names[i + 1], &names[high]);
    return (i + 1);
}

void quickSort(vector<string>& names, vector<double>& data, int low, int high)
{
    if (low < high)
    {
        int pi = partition(names, data, low, high);
        quickSort(names, data, low, pi - 1);
        quickSort(names, data, pi + 1, high);
    }
}

static int getMinSequenceLength(vector<string> names, vector<string> seqs)
{
    int currLen;
    int minLen = seqs[0].length();
    int minIndex = 0;
    for (int i = 1; i < seqs.size(); i++)
    {
        currLen = seqs[i].length();
        if (currLen < minLen)
        {
            minLen = currLen;
            minIndex = i;
        }
    }
    cout << "The shortest sequence is " << names[minIndex] << " - " << seqs[minIndex] << " and has a length of " << minLen << endl;
    return minLen;
}

static int checkFastaSameLength(vector<string> seqs)
{
    int length = seqs[0].length();
    for (int i = 1; i < seqs.size(); i++)
    {
        if (seqs[i].length() != length)
            return -1;
    }
    return length;
}

static void removeDisallowed(string& seq, const string& allowed) 
{
    unordered_set<char> allowedSet(allowed.begin(), allowed.end());
    seq.erase
    (
        remove_if(seq.begin(), seq.end(), [&](const char c) 
        {
            return !allowedSet.count(c);
        }
    ), seq.end());
}

static void translate(string& seq) {
    transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

    bool isProt = false;
    for (char c : seq)
    {
        bool found = false;
        for (char d : "ACTGU")
        {
            if (c == d)
            {
                found = true;
                break;
            }
        }
        if (!found)
        {
            isProt = true;
            break;
        }
    }
    if (isProt)
        return;

    replace(seq.begin(), seq.end(), 'U', 'T');
    string aaSeq;
    map<string, char> codonTable = getCodonTable();
    for (int i = 0; i < seq.length() - 2; i += 3)
    {
        string codon = seq.substr(i, 3);
        char aa = codonTable[codon];
        aaSeq.push_back(aa);
    }
    seq = aaSeq;
}

static vector<string> AAC(const string& seq, const string& seqName, const string& allowed)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    map<char, double> count;
    for(char c : allowed)
    {
        count[c] = 0;
    }
    int l = seq.length();
    for(int i = 0; i < l; i++) 
    {
        count[seq[i]]++;
    }
    int counter = 1;
    for(char c : allowed)
    {
        count[c] = count[c] / l;
        stringstream stream;
        stream << fixed << setprecision(7) << count[c];
        encoded.push_back(stream.str());
        counter++;
    }
    return encoded;
}

static vector<string> EAAC(const string& seq, const string& seqName, const string& allowed, vector<string> keys, int window, int windows)
{
    vector<string> encoded;
    encoded.push_back(seqName); 
    map<string, double> count;
    for(int i = 0; i < keys.size(); i++)
    {
        count[keys[i]] = 0;
    }
    for(int i = 0; i < windows; i++)
    {
        for(int j = 0; j < window; j++)
        {
            string key = to_string(i);
            key.push_back('-');
            key.push_back(seq[i + j]);
            count[key]++;
        }
    }
    int counter = 1;
    for(int i = 0; i < keys.size(); i++) 
    {
        count[keys[i]] = count[keys[i]] / window;
        stringstream stream;
        stream << fixed << setprecision(7) << count[keys[i]];
        encoded.push_back(stream.str());
        counter++;
    }
    return encoded;
}

static vector<string> CKSAAP(const string& seq, const string& seqName, const string& allowed, vector<string> keys, int gap)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    map<string, double> count;
    for(int i = 0; i < keys.size(); i++)
    {
        count[keys[i]] = 0;
    }
    int seqLen = seq.length();
    for(int i = 0; i < gap + 1; i++)
    {
        for(int j = 0; j < seqLen - i; j++)
        {
            string key = to_string(i);
            key.push_back('-');
            key.push_back(seq[j]);
            key.push_back(seq[i + j + 1]);
            count[key]++;
        }
    }
    int counter = 1;
    int gapCount = 0;
    int pairCount = 0;
    int allowedLen = allowed.length();
    int maxPairCount = allowedLen * allowedLen;
    for(int i = 0; i < keys.size(); i++) 
    {
        count[keys[i]] = count[keys[i]] / (seqLen - gapCount - 1);
        stringstream stream;
        stream << fixed << setprecision(7) << count[keys[i]];
        encoded.push_back(stream.str());
        counter++;
        pairCount++;
        if(pairCount == maxPairCount)
        {
            gapCount++;
            pairCount = 0;
        }
    }
    return encoded;
}

static vector<string> TPC(const string& seq, const string& seqName, const string& allowed, vector<string> keys)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    map<string, double> count;
    for(int i = 0; i < keys.size(); i++)
    {
        count[keys[i]] = 0;
    }
    int l = seq.length();
    for(int i = 0; i < l - 2; i++) 
    {
        string key(1, seq[i]);
        key.push_back(seq[i + 1]);
        key.push_back(seq[i + 2]);
        count[key]++;
    }
    int counter = 1;
    for(int i = 0; i < keys.size(); i++) 
    {
        count[keys[i]] = count[keys[i]] / (l - 2);
        stringstream stream;
        stream << fixed << setprecision(7) << count[keys[i]];
        encoded.push_back(stream.str());
        counter++;
    }
    return encoded;
}

static vector<string> DPC(const string& seq, const string& seqName, const string& allowed, vector<string> keys)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    map<string, double> count;
    for(int i = 0; i < keys.size(); i++)
    {
        count[keys[i]] = 0;
    }
    int l = seq.length();
    for(int i = 0; i < l - 1; i++) 
    {
        string key(1, seq[i]);
        key.push_back(seq[i + 1]);
        count[key]++;
    }
    int counter = 1;
    for(int i = 0; i < keys.size(); i++) 
    {
        count[keys[i]] = count[keys[i]] / (l - 1);
        stringstream stream;
        stream << fixed << setprecision(7) << count[keys[i]];
        encoded.push_back(stream.str());
        counter++;
    }
    return encoded;
}

static vector<string> DDE(const string& seq, const string& seqName, const string& allowed, vector<string> keys, map<string, double> tm, map<string, double> tvP)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    map<string, double> count;
    map<string, double> tv(tvP);
    int l = seq.length();
    for(int i = 0; i < keys.size(); i++)
    {
        count[keys[i]] = 0;
        tv[keys[i]] /= (l - 1);
    }

    for(int i = 0; i < l - 1; i++) 
    {
        string key(1, seq[i]);
        key.push_back(seq[i + 1]);
        count[key]++;
    }

    for(char c : allowed)
    {
        for(char d : allowed)
        {
            string key(1, c);
            key.push_back(d);
            count[key] /= (l - 1);
            count[key] = (count[key] - tm[key]) / sqrt(tv[key]);
        }
    }

    int counter = 1;
    for(int i = 0; i < keys.size(); i++) 
    {
        stringstream stream;
        stream << fixed << setprecision(7) << count[keys[i]];
        encoded.push_back(stream.str());
        counter++;
    }
    return encoded;
}

static vector<string> GAAC(const string& seq, const string& seqName, vector<string> keys, map<char, string> groups)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    map<string, double> count;
    for(int i = 0; i < keys.size(); i++)
    {
        count[keys[i]] = 0;
    }

    int l = seq.length();
    for(int i = 0; i < l; i++) 
    {
        count[groups[seq[i]]]++;
    }

    int counter = 1;
    for(int i = 0; i < keys.size(); i++) 
    {
        count[keys[i]] = count[keys[i]] / l;
        stringstream stream;
        stream << fixed << setprecision(7) << count[keys[i]];
        encoded.push_back(stream.str());
        counter++;
    }
    return encoded;
}

static vector<string> EGAAC(const string& seq, const string& seqName, vector<string> keys, map<char, string> groups, int window, int windows)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    map<string, double> count;
    for(int i = 0; i < keys.size(); i++)
    {
        count[keys[i]] = 0;
    }

    int l = seq.length();
    for(int i = 0; i < windows; i++)
    {
        for(int j = 0; j < window; j++) 
        {
            string key = to_string(i);
            key.push_back('-');
            key += groups[seq[i + j]];
            count[key]++;
        }
    }

    int counter = 1;
    for(int i = 0; i < keys.size(); i++) 
    {
        count[keys[i]] = count[keys[i]] / window;
        stringstream stream;
        stream << fixed << setprecision(7) << count[keys[i]];
        encoded.push_back(stream.str());
        counter++;
    }
    return encoded;
}

static vector<string> CKSAAGP(const string& seq, const string& seqName, vector<string> keys, int gap, map<char, string> groups, array<string, 5> groupStrings)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    map<string, double> count;
    for(int i = 0; i < keys.size(); i++)
    {
        count[keys[i]] = 0;
    }

    int seqLen = seq.length();
    for(int i = 0; i < gap + 1; i++)
    {
        for(int j = 0; j < seqLen - i; j++)
        {
            string key = to_string(i);
            key.push_back('-');
            key += groups[seq[j]];
            key.push_back('-');
            key += groups[seq[i + j + 1]];
            count[key]++;
        }
    }
    int counter = 1;
    int gapCount = 0;
    int pairCount = 0;
    int maxPairCount = 25;
    for(int i = 0; i < keys.size(); i++) 
    {
        count[keys[i]] = count[keys[i]] / (seqLen - gapCount - 1);
        stringstream stream;
        stream << fixed << setprecision(7) << count[keys[i]];
        encoded.push_back(stream.str());
        counter++;
        pairCount++;
        if(pairCount == maxPairCount)
        {
            gapCount++;
            pairCount = 0;
        }
    }
    return encoded;
}

static vector<string> GDPC(const string& seq, const string& seqName, vector<string> keys, map<char, string> groups, array<string, 5> groupStrings)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    map<string, double> count;
    for(int i = 0; i < keys.size(); i++)
    {
        count[keys[i]] = 0;
    }

    int l = seq.length();
    for(int i = 0; i < l - 1; i++) 
    {
        string key = groups[seq[i]];
        key.push_back('-');
        key += groups[seq[i + 1]];
        count[key]++;
    }
    int counter = 1;
    for(int i = 0; i < keys.size(); i++) 
    {
        count[keys[i]] = count[keys[i]] / (l - 1);
        stringstream stream;
        stream << fixed << setprecision(7) << count[keys[i]];
        encoded.push_back(stream.str());
        counter++;
    }
    return encoded;
}

static vector<string> GTPC(const string& seq, const string& seqName, vector<string> keys, map<char, string> groups, array<string, 5> groupStrings)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    map<string, double> count;
    for(int i = 0; i < keys.size(); i++)
    {
        count[keys[i]] = 0;
    }

    int l = seq.length();
    for(int i = 0; i < l - 2; i++) 
    {
        string key = groups[seq[i]];
        key.push_back('-');
        key += groups[seq[i + 1]];
        key.push_back('-');
        key += groups[seq[i + 2]];
        count[key]++;
    }
    int counter = 1;
    for(int i = 0; i < keys.size(); i++) 
    {
        count[keys[i]] = count[keys[i]] / (l - 2);
        stringstream stream;
        stream << fixed << setprecision(7) << count[keys[i]];
        encoded.push_back(stream.str());
        counter++;
    }
    return encoded;
}

static vector<string> binary(const string& seq, const string& seqName, vector<string> keys)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    map<string, double> count;
    for(int i = 0; i < keys.size(); i++)
    {
        count[keys[i]] = 0;
    }
    
    int counter = 0;
    for (char c : seq)
    {
        string key(1, c);
        key = to_string(counter) + key;
        count[key] = 1;
        counter++;
    }

    for(int i = 0; i < keys.size(); i++) 
    {
        stringstream stream;
        stream << fixed << setprecision(1) << count[keys[i]];
        encoded.push_back(stream.str());
    }
    return encoded;
}

static vector<string> Moran(const string& seq, const string& seqName, vector<string> keys, map<string, array<double, 20>> indices, vector<string> indexList, 
    int lag, map<string, double> means, map<char, int> order)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    map<string, double> average;
    int seqLen = seq.length();
    int count = 0;
    for (int in = 0; in < indexList.size(); in++)
    {
        string index = indexList[in];
        average[index] = 0;
        for (char c : seq)
        {
            average[index] += (indices[index][order[c]] / seqLen);
        }
        for (int i = 1; i <= lag; i++)
        {
            double num = 0;
            double den = 0;
            for (int j = 0; j < seqLen; j++)
            {
                char c1 = seq[j];
                char c2 = seq[j + i];
                double avg = average[index];
                if (j < seqLen - i)
                    num += (indices[index][order[c1]] - average[index]) * (indices[index][order[c2]] - average[index]);
                den += pow((indices[index][order[c1]] - average[index]), 2);
            }
            num /= (seqLen - i);
            den /= seqLen;
            double result = num / den;
            stringstream stream;
            stream << fixed << setprecision(7) << result;
            encoded.push_back(stream.str());
            count++;
        }
    }
    return encoded;
}

static vector<string> Geary(const string& seq, const string& seqName, vector<string> keys, map<string, array<double, 20>> indices, vector<string> indexList, 
    int lag, map<string, double> means, map<char, int> order)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    map<string, double> average;
    int seqLen = seq.length();
    int count = 0;
    for (int in = 0; in < indexList.size(); in++)
    {
        string index = indexList[in];
        average[index] = 0;
        for (char c : seq)
        {
            average[index] += (indices[index][order[c]] / seqLen);
        }
        for (int i = 1; i <= lag; i++)
        {
            double num = 0;
            double den = 0;
            for (int j = 0; j < seqLen; j++)
            {
                char c1 = seq[j];
                char c2 = seq[j + i];
                double avg = average[index];
                if (j < seqLen - i)
                    num += pow((indices[index][order[c1]] - indices[index][order[c2]]), 2);
                den += pow((indices[index][order[c1]] - average[index]), 2);
            }
            num /= (2*(seqLen - i));
            den /= (seqLen - 1);
            double result = num / den;
            stringstream stream;
            stream << fixed << setprecision(7) << result;
            encoded.push_back(stream.str());
            count++;
        }
    }
    return encoded;
}

static vector<string> NMB(const string& seq, const string& seqName, vector<string> keys, map<string, array<double, 20>> indices, vector<string> indexList, 
    int lag, map<string, double> means, map<char, int> order)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    int seqLen = seq.length();
    int count = 0;
    for (int in = 0; in < indexList.size(); in++)
    {
        string index = indexList[in];
        for (int i = 1; i <= lag; i++)
        {
            double num = 0;
            double den = 0;
            for (int j = 0; j < seqLen - i; j++)
            {
                char c1 = seq[j];
                char c2 = seq[j + i];
                num += (indices[index][order[c1]] * indices[index][order[c2]]);
            }
            den = seqLen - i;
            double result = num / den;
            stringstream stream;
            stream << fixed << setprecision(7) << result;
            encoded.push_back(stream.str());
            count++;
        }
    }
    return encoded;
}

static vector<string> CTDC(const string& seq, const string& seqName, const string& allowed, vector<string> keys, vector<map<string, string>> groups, 
    vector<string> properties)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    map<string, double> count;
    for(int i = 0; i < keys.size(); i++)
    {
        count[keys[i]] = 0;
    }
    for (int i = 0; i < properties.size(); i++)
    {
        for (int j = 0; j < groups.size(); j++)
        {
            string property = properties[i];
            for(char c : seq)
            {
                for (char d : groups[j][property])
                {
                    if (c == d)
                    {
                        string key = property;
                        key.push_back('-');
                        key += to_string(j + 1);
                        count[key]++;
                        goto endLoop;
                    }
                }
                endLoop:;
            }
        }
    }

    int l = seq.length();
    int counter = 1;
    for(int i = 0; i < keys.size(); i++) 
    {
        count[keys[i]] /= l;
        stringstream stream;
        stream << fixed << setprecision(7) << count[keys[i]];
        encoded.push_back(stream.str());
        counter++;
    }
    return encoded;
}

static vector<string> CTDT(const string& seq, const string& seqName, const string& allowed, vector<string> keys, vector<map<string, string>> groups, 
    vector<string> properties)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    map<string, double> count;
    for(int i = 0; i < keys.size(); i++)
    {
        count[keys[i]] = 0;
    }
    int seqLen = seq.length();
    for (int i = 0; i < properties.size(); i++)
    {
        for (int j = 0; j < groups.size(); j++)
        {
            for (int k = j + 1; k < groups.size(); k++)
            {
                string property = properties[i];
                for(int l = 0; l < seqLen - 1; l++)
                {
                    for (char c : groups[j][property])
                    {
                        for (char d : groups[k][property])
                        {
                            if ((seq[l] == c && seq[l + 1] == d) || (seq[l] == d && seq[l + 1] == c))
                            {
                                string key = property;
                                key.push_back('-');
                                key += to_string(j + 1);
                                key.push_back('-');
                                key += to_string(k + 1);
                                count[key]++;
                                goto endLoop;
                            }
                        }
                    }
                    endLoop:;
                }
            }
        }
    }

    int counter = 1;
    for(int i = 0; i < keys.size(); i++) 
    {
        count[keys[i]] /= (seqLen - 1);
        stringstream stream;
        stream << fixed << setprecision(7) << count[keys[i]];
        encoded.push_back(stream.str());
        counter++;
    }
    return encoded;
}

static vector<string> CTDD(const string& seq, const string& seqName, const string& allowed, vector<string> keys, vector<map<string, string>> groups, 
    vector<string> properties, array<int, 5> pcts)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    map<string, double> count;
    for(int i = 0; i < keys.size(); i++)
    {
        count[keys[i]] = 0;
    }
    int seqLen = seq.length();
    for (int i = 0; i < properties.size(); i++)
    {
        for (int j = 0; j < groups.size(); j++)
        {
            string property = properties[i];
            int counter = 0;
            for(char c : seq)
                for (char d : groups[j][property])
                    if (c == d)
                        counter++;
            double p25 = floor(0.25 * counter);
            double p50 = floor(0.50 * counter);
            double p75 = floor(0.75 * counter);
            double p100 = counter;
            map<int, double> cutoffs = {{0, 1}, {25, p25 > 1 ? p25 : 1}, {50, p50 > 1 ? p50 : 1}, {75, p75 > 1 ? p75 : 1}, {100, p100 > 1 ? p100 : 1}};
            for (auto &kvp : cutoffs)
            {
                string key = property;
                key.push_back('-');
                key += to_string(j + 1);
                key.push_back('-');
                key += to_string(kvp.first);
                counter = 0;
                for (int k = 0; k < seqLen; k++)
                {
                    for (char d : groups[j][property])
                    {
                        if (seq[k] == d)
                        {
                            counter++;
                            if (counter == kvp.second)
                            {
                                count[key] = (k + 1) * 100.0 / seqLen;
                                goto endLoop;
                            }
                        }
                    }
                }
                endLoop:;
                if (counter == 0)
                {
                    count[key] = 0;
                }
            }
        }
    }

    int counter = 1;
    for(int i = 0; i < keys.size(); i++) 
    {
        stringstream stream;
        stream << fixed << setprecision(7) << count[keys[i]];
        encoded.push_back(stream.str());
        counter++;
    }
    return encoded;
}

static vector<string> CT(const string& seq, const string& seqName, const string& allowed, vector<string> keys, array<string, 7> groups)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    map<string, double> count;
    for(int i = 0; i < keys.size(); i++)
    {
        count[keys[i]] = 0;
    }
    int l = seq.length();
    for (int i = 0; i < l - 2; i++) 
    {
        int g1 = 0;
        int g2 = 0;
        int g3 = 0;
        for (int j = 0; j < 7; j++)
        {
            for (char c : groups[j])
            {
                if (seq[i] == c)
                    g1 = j + 1;
                if (seq[i + 1] == c)
                    g2 = j + 1;
                if (seq[i + 2] == c)
                    g3 = j + 1;
            }
        }
        string key = to_string(g1);
        key.push_back('-');
        key += to_string(g2);
        key.push_back('-');
        key += to_string(g3);
        count[key]++;
    }
    double min = l;
    double max = 0;
    for (auto &kvp : count)
    {
        min = kvp.second < min ? kvp.second : min;
        max = kvp.second > max ? kvp.second : max;
    }
    int counter = 1;
    for(int i = 0; i < keys.size(); i++) 
    {
        count[keys[i]] = (count[keys[i]] - min) / max;
        stringstream stream;
        stream << fixed << setprecision(7) << count[keys[i]];
        encoded.push_back(stream.str());
        counter++;
    }
    return encoded;
}

static vector<string> KSCT(const string& seq, const string& seqName, const string& allowed, vector<string> keys, array<string, 7> groups, int ks)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    map<string, double> count;
    for(int i = 0; i < keys.size(); i++)
    {
        count[keys[i]] = 0;
    }
    int l = seq.length();
    for (int currK = 0; currK <= ks; currK++)
    {
        for (int i = 0; i < l - 2; i++) 
        {
            int g1 = 0;
            int g2 = 0;
            int g3 = 0;
            for (int j = 0; j < 7; j++)
            {
                for (char c : groups[j])
                {
                    if (seq[i] == c)
                        g1 = j + 1;
                    if (seq[i + (1 * (currK + 1))] == c)
                        g2 = j + 1;
                    if (seq[i + (2 * (currK + 1))] == c)
                        g3 = j + 1;
                }
            }
            string key = to_string(g1);
            key.push_back('-');
            key += to_string(g2);
            key.push_back('-');
            key += to_string(g3);
            key.push_back('-');
            key += to_string(currK);
            count[key]++;
        }
        double min = l;
        double max = 0;
        for (int i = currK * 343; i < (currK + 1) * 343; i++)
        {
            min = count[keys[i]] < min ? count[keys[i]] : min;
            max = count[keys[i]] > max ? count[keys[i]] : max;
        }
        for (int i = currK * 343; i < (currK + 1) * 343; i++)
        {
            count[keys[i]] = (count[keys[i]] - min) / max;
        }
    } 
    
    int counter = 1;
    for(int i = 0; i < keys.size(); i++) 
    {
        stringstream stream;
        stream << fixed << setprecision(7) << count[keys[i]];
        encoded.push_back(stream.str());
        counter++;
    }
    return encoded;
}

static vector<string> SOCNumber(const string& seq, const string& seqName, const string& allowed, vector<string> keys, map<char, int> order, 
    map<char, array<double, 20>> swData, map<char, array<double, 20>> gData, int lag)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    map<string, double> count;
    for (int i = 0; i < keys.size(); i++)
    {
        count[keys[i]] = 0;
    }
    int l = seq.length();
    for (int i = 1; i <= lag; i++)
    {
        double sw = 0;
        double g = 0;
        for (int j = 0; j < l - i; j++)
        {
            sw += pow(swData[seq[j]][order[seq[j + i]]], 2);
            g += pow(gData[seq[j]][order[seq[j + i]]], 2);
        }
        sw /= (l - i);
        g /= (l - i);

        string key(1, '-');
        key += to_string(i);
        count["sw" + key] = sw;
        count["g" + key] = g;
    }

    int counter = 1;
    for (int i = 0; i < keys.size(); i++) 
    {
        stringstream stream;
        stream << fixed << setprecision(7) << count[keys[i]];
        encoded.push_back(stream.str());
        counter++;
    }
    return encoded;
}

static vector<string> QSOrder(const string& seq, const string& seqName, const string& allowed, vector<string> keys, map<char, int> order, 
    map<char, array<double, 20>> swData, map<char, array<double, 20>> gData, int lag, double weight)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    map<string, double> count;
    for (int i = 0; i < keys.size(); i++)
    {
        count[keys[i]] = 0;
    }
    int l = seq.length();
    for (int i = 1; i <= lag; i++)
    {
        double sw = 0;
        double g = 0;
        for (int j = 0; j < l - i; j++)
        {
            sw += pow(swData[seq[j]][order[seq[j + i]]], 2);
            g += pow(gData[seq[j]][order[seq[j + i]]], 2);
        }

        string key(1, '-');
        key += to_string(i);
        count["sw" + key] = sw;
        count["g" + key] = g;
    }

    double sumSW = 0;
    double sumG = 0;
    for(int i = 1; i <= lag; i++)
    {
        sumSW += count["sw-" + to_string(i)];
        sumG += count["g-" + to_string(i)];
    }

    map<char, int> counts;
    for (char c : allowed)
    {
        counts[c] = 0;
    }
    for (char c : seq)
    {
        counts[c]++;
    }
    for (char c : allowed)
    {
        string key(1, '-');
        key.push_back(c);
        count["sw" + key] = counts[c] / (1 + weight * sumSW);
        count["g" + key] = counts[c] / (1 + weight * sumG);
    }
    for (int i = 1; i <= lag; i++)
    {
        string key(1, '-');
        key += to_string(i);
        count["sw" + key] *= weight / (1 + weight * sumSW);
        count["g" + key] *= weight / (1 + weight * sumG);
    }

    int counter = 1;
    for (int i = 0; i < keys.size(); i++) 
    {
        stringstream stream;
        stream << fixed << setprecision(7) << count[keys[i]];
        encoded.push_back(stream.str());
        counter++;
    }
    return encoded;
}

static vector<string> PAAC(const string& seq, const string& seqName, const string& allowed, vector<string> keys, int lag, map<char, int> order,
    map<string, array<double, 20>> paac, double weight)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    map<string, double> count;
    for (int i = 0; i < keys.size(); i++)
    {
        count[keys[i]] = 0;
    }
    int l = seq.length();
    double globalSum = 0;
    for (int i = 1; i <= lag; i++)
    {
        double res = 0;
        for (int j = 0; j < l - i; j++)
        {
            for (auto &kvp : paac)
            {
                res += pow(kvp.second[order[seq[j]]] - kvp.second[order[seq[j + i]]], 2) / 3; // 3 for H1, H2 and SCM
            }
        }
        res /= (l - i);
        globalSum += res;
        string key = to_string(i);
        count[key] = res;
    }
    map<char, int> counts;
    for (char c : allowed)
    {
        counts[c] = 0;
    }
    for (char c : seq)
    {
        counts[c]++;
    }

    for (char c : allowed)
    {
        string key(1, c);
        count[key] = counts[c] / (1 + weight * globalSum);
    }
    for (int i = 1; i <= lag; i++)
    {
        string key = to_string(i);
        count[key] = (weight * count[key]) / (1 + weight * globalSum);
    }

    int counter = 1;
    for (int i = 0; i < keys.size(); i++) 
    {
        stringstream stream;
        stream << fixed << setprecision(7) << count[keys[i]];
        encoded.push_back(stream.str());
        counter++;
    }
    return encoded;
}

static vector<string> APAAC(const string& seq, const string& seqName, const string& allowed, vector<string> keys, int lag, map<char, int> order,
    map<string, array<double, 20>> paac, string properties[], double weight)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    map<string, double> count;
    for (int i = 0; i < keys.size(); i++)
    {
        count[keys[i]] = 0;
    }
    
    int l = seq.length();
    double globalSum = 0;
    for (int i = 1; i <= lag; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            string property = properties[j];
            double res = 0;
            for (int k = 0; k < l - i; k++)
                res += paac[property][order[seq[k]]] * paac[property][order[seq[k + i]]];
            res /= (l - i);
            globalSum += res;
            string key = property;
            key.push_back('-');
            key += to_string(i);
            count[key] = res;
        }
    }

    map<char, int> counts;  
    for (char c : allowed)
    {
        counts[c] = 0;
    }
    for (char c : seq)
    {
        counts[c]++;
    }

    for (char c : allowed)
    {
        string key(1, c);
        count[key] = counts[c] / (1 + weight * globalSum);
    }
    for (int i = 1; i <= lag; i++)
    {
        for (auto &kvp : paac)
        {
            string key = kvp.first;
            key.push_back('-');
            key += to_string(i);
            count[key] = (weight * count[key]) / (1 + weight * globalSum);
        }
    }

    int counter = 1;
    for (int i = 0; i < keys.size(); i++) 
    {
        stringstream stream;
        stream << fixed << setprecision(7) << count[keys[i]];
        encoded.push_back(stream.str());
        counter++;
    }
    return encoded;
}

static vector<string> AAINDEX(const string& seq, const string& seqName, vector<string> keys, map<char, int> order, map<string, 
    array<double, 20>> indexData, vector<string> indexList)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    map<string, double> count;
    for (int i = 0; i < keys.size(); i++)
    {
        count[keys[i]] = 0;
    }
    
    int l = seq.length();
    int counter = 1;
    for (int i = 0; i < l; i++)
    {
        for (string index : indexList)
        {
            stringstream stream;
            stream << fixed << setprecision(7) << indexData[index][order[seq[i]]];
            encoded.push_back(stream.str());
            counter++;
        }
    }
    return encoded;
}

static vector<string> BLOSUM62(const string& seq, const string& seqName, vector<string> keys, map<char, array<double, 20>> blosum)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    int l = seq.length();
    int counter = 1;
    for (int i = 0; i < l; i++)
    {
        for (int j = 0; j < 20; j++)
        {
            stringstream stream;
            stream << fixed << setprecision(7) << blosum[seq[i]][j];
            encoded.push_back(stream.str());
            counter++;
        }
    }
    return encoded;
}

static vector<string> ZSCALE(const string& seq, const string& seqName, vector<string> keys, map<char, array<double, 5>> zscale)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    int l = seq.length();
    int counter = 1;
    for (int i = 0; i < l; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            stringstream stream;
            stream << fixed << setprecision(7) << zscale[seq[i]][j];
            encoded.push_back(stream.str());
            counter++;
        }
    }
    return encoded;
}

static vector<string> SSEB(const string& seq, const string& seqName, const string& allowed, vector<string> keys, ifstream &file, map<char, array<int, 3>> types, const string& type)
{
    vector<string> encoded;
    encoded.push_back(seqName);

    string aas;
    string elements;
    // Read PSSM file
    string line;
    getline(file, line); // Skip the first line
    if (type.compare("ss2") == 0)
        getline(file, line); // Skip the second line if .ss2

    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for(string line; iss >> line; counter++) 
        {
            if (counter == 1)
            {
                bool exists = false;
                for (char c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
                aas.push_back(line[0]);
            }
            else if (counter == 2)
            {
                elements.push_back(line[0]);
                break;
            };
        }
        endLoop:;
    }
    size_t found = aas.find(seq);
    if (found == string::npos)
        cout << "Error: Sequence " << seqName << " not found in the file.";
    else
    {
        int l = seq.length();
        for (int i = found; i < found + l; i++)
        {
            char c = elements[i];
            if (types.find(c) != types.end())
            {
                array<int, 3> vals = types[c];
                for (int val : vals)
                {
                    stringstream stream;
                    stream << fixed << setprecision(7) << (val * 1.0);
                    encoded.push_back(stream.str());
                } 
            }
            else 
            {
                encoded.push_back("0.0000000");
                encoded.push_back("0.0000000");
                encoded.push_back("0.0000000");
            }
        }
    }
                
    return encoded;
}

static vector<string> SSEC(const string& seq, const string& seqName, const string& allowed, vector<string> keys, ifstream &file, array<char, 3> types, const string& type)
{
    vector<string> encoded;
    encoded.push_back(seqName);

    string aas;
    string elements;
    // Read PSSM file
    string line;
    getline(file, line); // Skip the first line
    if (type.compare("ss2") == 0)
        getline(file, line); // Skip the second line if .ss2
    
    map<char, int> typesCount;
    for (char c : types)
    {
        typesCount[c] = 0;
    }
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for(string line; iss >> line; counter++) 
        {
            if (counter == 1)
            {
                bool exists = false;
                for (char c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
                aas.push_back(line[0]);
            }
            else if (counter == 2)
            {
                elements.push_back(line[0]);
                break;
            };
        }
        endLoop:;
    }

    size_t found = aas.find(seq);
    if (found == string::npos)
        cout << "Error: Sequence " << seqName << " not found in the file.";
    else
    {
        int l = seq.length();
        for (int i = found; i < found + l; i++)
        {
            char c = elements[i];
            if (typesCount.find(c) != typesCount.end())
                typesCount[c]++;
        }
        for (char c : types)
        {
            double val = typesCount[c] / (l * 1.0);
            stringstream stream;
            stream << fixed << setprecision(7) << val;
            encoded.push_back(stream.str());
        } 
    }
    
    return encoded;
}

static vector<string> SSPB(const string& seq, const string& seqName, const string& allowed, vector<string> keys, int n, ifstream &file, array<char, 3> types, const string& type)
{
    vector<string> encoded;
    encoded.push_back(seqName);

    string aas;
    map<char, vector<double>> probs;
    // Read PSSM file
    string line;
    getline(file, line); // Skip the first line
    if (type.compare("ss2") == 0)
        getline(file, line); // Skip the second line if .ss2
    
    map<string, double> count;
    for (string key : keys)
    {
        count[key] = 0;
    }
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for(string line; iss >> line; counter++) 
        {
            char currType;
            if (counter == 1)
            {
                bool exists = false;
                for (char c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
                aas.push_back(line[0]);
            }
            else if ((counter == 3 && type.compare("ss2") == 0) || (counter == 7 && type.compare("spXout") == 0))
            {
                double val = stof(line);
                probs['H'].push_back(val);
            }
            else if ((counter == 4 && type.compare("ss2") == 0) || (counter == 5 && type.compare("spXout") == 0))
            {
                double val = stof(line);
                probs['E'].push_back(val);
            }
            else if ((counter == 5 && type.compare("ss2") == 0) || (counter == 6 && type.compare("spXout") == 0))
            {
                double val = stof(line);
                probs['C'].push_back(val);
                break;
            }
        }
        endLoop:;
    }

    size_t found = aas.find(seq);
    if (found == string::npos)
        cout << "Error: Sequence " << seqName << " not found in the file.";
    else
    {
        int l = seq.length();
        for (char c : types)
        {
            for (char d : types)
            {
                string key(1, c);
                key.push_back(d);
                for (int i = found; i < found + l - n; i++)
                {
                    count[key] += probs[c][i] * probs[d][i + n];
                }
            }
        }

        for (string key : keys)
        {
            double val = count[key];
            stringstream stream;
            stream << fixed << setprecision(7) << val;
            encoded.push_back(stream.str());
        } 
    }
    
    return encoded;
}

static vector<string> SSPAC(const string& seq, const string& seqName, const string& allowed, vector<string> keys, int n, ifstream &file, array<char, 3> types, const string& type)
{
    vector<string> encoded;
    encoded.push_back(seqName);

    string aas;
    map<char, vector<double>> probs;
    // Read PSSM file
    string line;
    getline(file, line); // Skip the first line
    if (type.compare("ss2") == 0)
        getline(file, line); // Skip the second line if .ss2
    
    map<string, double> count;
    for (string key : keys)
    {
        count[key] = 0;
    }
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for(string line; iss >> line; counter++) 
        {
            char currType;
            if (counter == 1)
            {
                bool exists = false;
                for (char c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
                aas.push_back(line[0]);
            }
            else if ((counter == 3 && type.compare("ss2") == 0) || (counter == 7 && type.compare("spXout") == 0))
            {
                double val = stof(line);
                probs['H'].push_back(val);
            }
            else if ((counter == 4 && type.compare("ss2") == 0) || (counter == 5 && type.compare("spXout") == 0))
            {
                double val = stof(line);
                probs['E'].push_back(val);
            }
            else if ((counter == 5 && type.compare("ss2") == 0) || (counter == 6 && type.compare("spXout") == 0))
            {
                double val = stof(line);
                probs['C'].push_back(val);
                break;
            }
        }
        endLoop:;
    }

    size_t found = aas.find(seq);
    if (found == string::npos)
        cout << "Error: Sequence " << seqName << " not found in the file.";
    else
    {
        int l = seq.length();
        for (int i = 1; i <= n; i++)
        {
            string iString = to_string(i);
            for (char c : types)
            {
                string key = iString;
                key.push_back(c);
                for (int j = found; j < found + l - i; j++)
                { 
                    count[key] += probs[c][j] * probs[c][j + i];
                }
            }
        }

        for (string key : keys)
        {
            double val = count[key] / l;
            stringstream stream;
            stream << fixed << setprecision(7) << val;
            encoded.push_back(stream.str());
        } 
    }
    
    return encoded;
}

static vector<string> Disorder(const string& seq, const string& seqName, const string& allowed, vector<string> keys, ifstream &file)
{
    vector<string> encoded;
    encoded.push_back(seqName);

    // Read disorder file
    string line;
    string aas;
    vector<double> scores;
    getline(file, line);
    while (line[0] != '-')
    {
        getline(file, line); 
    }
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for(string line; iss >> line; counter++) 
        {
            if (line[0] == '=')
                break; 
            else if (counter == 1)
            {
                bool exists = false;
                for (char c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
                aas.push_back(line[0]);
            }
            else if (counter == 2)
            {
                double val = stof(line);
                scores.push_back(val);
                break;
            };
        }
        endLoop:;
    }
    size_t found = aas.find(seq);
    if (found == string::npos)
        cout << "Error: Sequence " << seqName << " not found in the file.";
    else
    {
        int l = seq.length();
        for (int i = found; i < found + l; i++)
        {
            double val = scores[i];
            stringstream stream;
            stream << fixed << setprecision(7) << val;
            encoded.push_back(stream.str());
        }
    }

    return encoded;
}

static vector<string> DisorderC(const string& seq, const string& seqName, const string& allowed, vector<string> keys, ifstream &file)
{
    vector<string> encoded;
    encoded.push_back(seqName);

    // Read disorder file
    string line;
    string aas;
    string types;
    getline(file, line);
    while (line[0] != '-')
    {
        getline(file, line); 
    }
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for(string line; iss >> line; counter++) 
        {
            if (line[0] == '=')
                break; 
            else if (counter == 1)
            {
                bool exists = false;
                for (char c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
                aas.push_back(line[0]);
            }
            else if (counter == 3)
            {
                char val = line[0];
                types.push_back(val);
                break;
            };
        }
        endLoop:;
    }
    size_t found = aas.find(seq);
    if (found == string::npos)
        cout << "Error: Sequence " << seqName << " not found in the file.";
    else
    {
        int l = seq.length();
        double order = 0;
        double disorder = 0;
        for (int i = found; i < found + l; i++)
        {
            if (types[i] == 'D')
                disorder++;
            else
                order++;
        }
        disorder /= l;
        order /= l;
        stringstream stream;
        stream << fixed << setprecision(7) << disorder;
        encoded.push_back(stream.str());
        stringstream stream2;
        stream2 << fixed << setprecision(7) << order;
        encoded.push_back(stream2.str());
    }

    return encoded;
}

static vector<string> DisorderB(const string& seq, const string& seqName, const string& allowed, vector<string> keys, ifstream &file)
{
    vector<string> encoded;
    encoded.push_back(seqName);

    // Read disorder file
    string line;
    string aas;
    string types;
    getline(file, line);
    while (line[0] != '-')
    {
        getline(file, line); 
    }
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for(string line; iss >> line; counter++) 
        {
            if (line[0] == '=')
                break; 
            else if (counter == 1)
            {
                bool exists = false;
                for (char c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
                aas.push_back(line[0]);
            }
            else if (counter == 3)
            {
                char val = line[0];
                types.push_back(val);
                break;
            };
        }
        endLoop:;
    }
    size_t found = aas.find(seq);
    if (found == string::npos)
        cout << "Error: Sequence " << seqName << " not found in the file.";
    else
    {
        int l = seq.length();
        for (int i = found; i < found + l; i++)
        {
            if (types[i] == 'D')
            {
                encoded.push_back("0.0000000");
                encoded.push_back("1.0000000");
            }
            else
            {
                encoded.push_back("1.0000000");
                encoded.push_back("0.0000000");
            }
        }
    }

    return encoded;
}

static vector<string> TA(const string& seq, const string& seqName, const string& allowed, vector<string> keys, ifstream &file)
{
    vector<string> encoded;
    encoded.push_back(seqName);

    // Read disorder file
    string line;
    string aas;
    vector<double> phis;
    vector<double> psis;
    getline(file, line); // Skip first line
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for(string line; iss >> line; counter++) 
        {
            if (counter == 1)
            {
                bool exists = false;
                for (char c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
                aas.push_back(line[0]);
            }
            else if (counter == 3)
            {
                double val = stof(line);
                phis.push_back(val);
            }
            else if (counter == 4)
            {
                double val = stof(line);
                psis.push_back(val);
                break;
            }
        }
        endLoop:;
    }
    size_t found = aas.find(seq);
    if (found == string::npos)
        cout << "Error: Sequence " << seqName << " not found in the file.";
    else
    {
        int l = seq.length();
        for (int i = found; i < found + l; i++)
        {
            stringstream stream;
            stream << fixed << setprecision(7) << phis[i];
            encoded.push_back(stream.str());
            stringstream stream2;
            stream2 << fixed << setprecision(7) << psis[i];
            encoded.push_back(stream2.str());
        }
    }

    return encoded;
}

static vector<string> TAC(const string& seq, const string& seqName, const string& allowed, vector<string> keys, ifstream &file)
{
    vector<string> encoded;
    encoded.push_back(seqName);

    double pi = 3.14159265359;
    
    // Read disorder file
    string line;
    string aas;
    vector<double> phiSin;
    vector<double> phiCos;
    vector<double> psiSin;
    vector<double> psiCos;
    getline(file, line); // Skip first line
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for(string line; iss >> line; counter++) 
        {
            if (counter == 1)
            {
                bool exists = false;
                for (char c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
                aas.push_back(line[0]);
            }
            else if (counter == 3)
            {
                double val = stof(line) * (pi / 180);
                phiSin.push_back(sin(val));
                phiCos.push_back(cos(val));
            }
            else if (counter == 4)
            {
                double val = stof(line) * (pi / 180);
                psiSin.push_back(sin(val));
                psiCos.push_back(cos(val));
                break;
            }
        }
        endLoop:;
    }
    size_t found = aas.find(seq);
    if (found == string::npos)
        cout << "Error: Sequence " << seqName << " not found in the file.";
    else
    {
        int l = seq.length();
        double phiSinVal = 0;
        double phiCosVal = 0;
        double psiSinVal = 0;
        double psiCosVal = 0;
        for (int i = found; i < found + l; i++)
        {
            phiSinVal += phiSin[i];
            phiCosVal += phiCos[i];
            psiSinVal += psiSin[i];
            psiCosVal += psiCos[i];
        }
        stringstream stream;
        stream << fixed << setprecision(7) << phiSinVal / l;
        encoded.push_back(stream.str());
        stringstream stream2;
        stream2 << fixed << setprecision(7) << phiCosVal / l;
        encoded.push_back(stream2.str());
        stringstream stream3;
        stream3 << fixed << setprecision(7) << psiSinVal / l;
        encoded.push_back(stream3.str());
        stringstream stream4;
        stream4 << fixed << setprecision(7) << psiCosVal / l;
        encoded.push_back(stream4.str());
    }

    return encoded;
}

static vector<string> TAB(const string& seq, const string& seqName, const string& allowed, vector<string> keys, int n, ifstream &file)
{
    vector<string> encoded;
    encoded.push_back(seqName);

    double pi = 3.14159265359;
    
    // Read disorder file
    string line;
    string aas;
    vector<double> phiSin;
    vector<double> phiCos;
    vector<double> psiSin;
    vector<double> psiCos;
    getline(file, line); // Skip first line
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for(string line; iss >> line; counter++) 
        {
            if (counter == 1)
            {
                bool exists = false;
                for (char c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
                aas.push_back(line[0]);
            }
            else if (counter == 3)
            {
                double val = stof(line) * (pi / 180);
                phiSin.push_back(sin(val));
                phiCos.push_back(cos(val));
            }
            else if (counter == 4)
            {
                double val = stof(line) * (pi / 180);
                psiSin.push_back(sin(val));
                psiCos.push_back(cos(val));
                break;
            }
        }
        endLoop:;
    }
    size_t found = aas.find(seq);
    if (found == string::npos)
        cout << "Error: Sequence " << seqName << " not found in the file.";
    else
    {
        map<string, double> values;
        int l = seq.length();
        for (int i = found; i < found + l - n; i++)
        {
            values["phiSin-phiSin"] = phiSin[i] * phiSin[i + n];
            values["phiSin-phiCos"] = phiSin[i] * phiCos[i + n];
            values["phiSin-psiSin"] = phiSin[i] * psiSin[i + n];
            values["phiSin-psiCos"] = phiSin[i] * psiCos[i + n];
            values["phiCos-phiCos"] = phiCos[i] * phiCos[i + n];
            values["phiCos-psiSin"] = phiCos[i] * psiSin[i + n];
            values["phiCos-psiCos"] = phiCos[i] * psiCos[i + n];
            values["psiSin-psiSin"] = psiSin[i] * psiSin[i + n];
            values["psiSin-psiCos"] = psiSin[i] * psiCos[i + n];
            values["psiCos-psiCos"] = psiCos[i] * psiCos[i + n];
        }

        for (string key : keys)
        {
            values[key] /= l;
            stringstream stream;
            stream << fixed << setprecision(7) << values[key];
            encoded.push_back(stream.str());
        }
    }

    return encoded;
}

static vector<string> TAAC(const string& seq, const string& seqName, const string& allowed, vector<string> keys, int n, ifstream &file)
{
    vector<string> encoded;
    encoded.push_back(seqName);

    double pi = 3.14159265359;
    
    // Read disorder file
    string line;
    string aas;
    vector<double> phiSin;
    vector<double> phiCos;
    vector<double> psiSin;
    vector<double> psiCos;
    getline(file, line); // Skip first line
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for(string line; iss >> line; counter++) 
        {
            if (counter == 1)
            {
                bool exists = false;
                for (char c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
                aas.push_back(line[0]);
            }
            else if (counter == 3)
            {
                double val = stof(line) * (pi / 180);
                phiSin.push_back(sin(val));
                phiCos.push_back(cos(val));
            }
            else if (counter == 4)
            {
                double val = stof(line) * (pi / 180);
                psiSin.push_back(sin(val));
                psiCos.push_back(cos(val));
                break;
            }
        }
        endLoop:;
    }
    size_t found = aas.find(seq);
    if (found == string::npos)
        cout << "Error: Sequence " << seqName << " not found in the file.";
    else
    {
        map<string, double> values;
        int l = seq.length();
        for (int i = 1; i <= n; i++)
        {
            string iString = to_string(i);
            for (int j = found; j < found + l - i; j++)
            {
                values[iString + "-phiSin"] = phiSin[i] * phiSin[i + n];
                values[iString + "-phiCos"] = phiCos[i] * phiCos[i + n];
                values[iString + "-psiSin"] = psiSin[i] * psiSin[i + n];
                values[iString + "-psiCos"] = psiCos[i] * psiCos[i + n];
            }
        }

        for (string key : keys)
        {
            values[key] /= l;
            stringstream stream;
            stream << fixed << setprecision(7) << values[key];
            encoded.push_back(stream.str());
        }
    }

    return encoded;
}

static vector<string> ASA(const string& seq, const string& seqName, const string& allowed, vector<string> keys, ifstream &file)
{
    vector<string> encoded;
    encoded.push_back(seqName);

    // Read disorder file
    string line;
    string aas;
    vector<double> asas;
    getline(file, line); // Skip first line
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for(string line; iss >> line; counter++) 
        {
            if (counter == 1)
            {
                bool exists = false;
                for (char c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
                aas.push_back(line[0]);
            }
            else if (counter == 10)
            {
                double val = stof(line);
                asas.push_back(val);
                break;
            }
        }
        endLoop:;
    }
    size_t found = aas.find(seq);
    if (found == string::npos)
        cout << "Error: Sequence " << seqName << " not found in the file.";
    else
    {
        int l = seq.length();
        for (int i = found; i < found + l; i++)
        {
            asas[i] /= l;
            stringstream stream;
            stream << fixed << setprecision(7) << asas[i];
            encoded.push_back(stream.str());
        }
    }

    return encoded;
}

static vector<string> KNNpeptide(const string& seq, const string& seqName, vector<string> keys, map<char, array<double, 20>> blosum, map<char, int> order,
    vector<string> trainNames, vector<string> trainSeqs, map<string, string> labelsData, set<string> labelsP, vector<int> kNums)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    
    vector<string> labels;
    vector<double> distances;
    for (int i = 0; i < trainSeqs.size(); i++)
    {
        if (seqName != trainNames[i])
        {
            labels.push_back(labelsData[trainNames[i]]);
            double minValue = -4.0;
            double maxValue = 11.0;
            double similarity = 0;
            for (int j = 0; j < trainSeqs[i].length(); j++)
            {
                similarity += (blosum[trainSeqs[i][j]][order[seq[j]]] - minValue) / (maxValue - minValue);
            }
            similarity = 1 - similarity / trainSeqs[i].length();
            distances.push_back(similarity);
        }
    }

    quickSort(labels, distances, 0, distances.size() - 1);
    for (int num : kNums)
    {
        map<string, double> contents;
        for (string label : labelsP)
        {
            contents[label] = 0;
        }
        for (int i = 0; i < num; i++)
        {
            contents[labels[i]]++;
        }
        for (string label : labelsP)
        {
            double val = contents[label] /= num;
            stringstream stream;
            stream << fixed << setprecision(7) << val;
            encoded.push_back(stream.str());
        }
    }
    
    return encoded;
}

static vector<string> KNNprotein(const string& seq, const string& seqName, vector<string> keys, map<char, array<double, 20>> blosum, map<char, int> order,
    vector<string> trainNames, vector<string> trainSeqs, map<string, string> labelsData, set<string> labelsP, vector<int> kNums)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    
    vector<string> labels;
    vector<double> distances;
    for (int i = 0; i < trainSeqs.size(); i++)
    {
        if (seqName != trainNames[i])
        {
            labels.push_back(labelsData[trainNames[i]]);
            double similarity = 0;
            array<string, 2> nw = calculateNeedlemanWunsch(trainSeqs[i], seq, blosum, order, -10, -1);
            for (int j = 0; j < nw[0].length(); j++)
            {
                if (nw[0][j] == nw[1][j])
                    similarity++;
            }
            similarity *= 2 / (trainSeqs[i].length() + seq.length());
            distances.push_back(similarity);
        }
    }

    quickSort(labels, distances, 0, distances.size() - 1);
    for (int num : kNums)
    {
        map<string, double> contents;
        for (string label : labelsP)
        {
            contents[label] = 0;
        }
        for (int i = 0; i < num; i++)
        {
            contents[labels[i]]++;
        }
        for (string label : labelsP)
        {
            double val = contents[label] /= num;
            stringstream stream;
            stream << fixed << setprecision(7) << val;
            encoded.push_back(stream.str());
        }
    }

    return encoded;
}

static vector<string> AAPAS(const string& seq, const string& seqName, const string& allowed, vector<string> keys, map<string, double> freqs)
{
    vector<string> encoded;
    encoded.push_back(seqName);

    map<string, double> count;
    
    int l = seq.length();
    for(int i = 0; i < l - 1; i++) 
    {
        string key(1, seq[i]);
        key.push_back(seq[i + 1]);
        count[key] += freqs[key];
    }
    for(int i = 0; i < keys.size(); i++) 
    {
        stringstream stream;
        stream << fixed << setprecision(7) << count[keys[i]];
        encoded.push_back(stream.str());
    }
    return encoded;
}

static vector<string> TVD(const string& seq, const string& seqName, const string& allowed, vector<string> keys, map<char, array<double, 10>> tvd)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    
    int l = seq.length();
    for(int i = 0; i < l; i++) 
    {
        for(int j = 0; j < 10; j++)
        {
            stringstream stream;
            stream << fixed << setprecision(7) << tvd[seq[i]][j];
            encoded.push_back(stream.str());
        }
    }
    return encoded;
}

static vector<string> CMV(const string& seq, const string& seqName, const string& allowed, vector<string> keys)
{
    vector<string> encoded;
    encoded.push_back(seqName);

    map<char, double> count;
    for (char c : allowed)
    {
        count[c] = 0;
    }

    int l = seq.length();
    for (int i = 0; i < l; i++)
    {
        count[seq[i]] += i + 1;
    }

    for(int i = 0; i < l; i++) 
    {
        for(char c : allowed)
        {
            count[c] /= (l * (l - 1));
            stringstream stream;
            stream << fixed << setprecision(7) << count[c];
            encoded.push_back(stream.str());
        }
    }
    return encoded;
}

static vector<string> EBGW(const string& seq, const string& seqName, const string& allowed, vector<string> keys, map<char, array<int, 3>> groups, int n)
{
    vector<string> encoded;
    encoded.push_back(seqName);

    array<vector<double>, 3> count;

    int l = seq.length();

    for (int k = 1; k <= n; k++)
    {
        for(int i = 0; i < 3; i++)
        {
            count[i].push_back(0);
        }
        int subSize = floor(k * l / (n * 1.0));
        if (subSize == 0)
        {
            for(int i = 0; i < 3; i++)
            {
                count[i][k - 1] = 0;
            }
        }
        else
        {
            double subSizeD = subSize * 1.0;
            string sub = seq.substr(0, subSize);
            for (char c : sub)
            {
                for(int i = 0; i < 3; i++)
                {
                    count[i][k - 1] += groups[c][i];
                }
            }
            for(int i = 0; i < 3; i++)
            {
                count[i][k - 1] /= subSizeD;
            }
        }
    }

    for(vector<double> v : count)
    {
        for (double d : v)
        {
            stringstream stream;
            stream << fixed << setprecision(7) << d;
            encoded.push_back(stream.str());
        }
    }

    return encoded;
}

static vector<string> PSSM(const string& seq, const string& seqName, const string& allowed, vector<string> keys, ifstream &file)
{
    vector<string> encoded;
    encoded.push_back(seqName);

    // Read PSSM file
    string line;
    getline(file, line); // Skip the first line
    getline(file, line); // Skip the second line
    getline(file, line); // Skip the third line
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for(string line; iss >> line && counter < 22;) 
        {
            if(counter == 1)
            {
                bool exists = false;
                for (char c : allowed)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
            }
            if (counter >= 2)
            {
                double val = stof(line);
                stringstream stream;
                stream << fixed << setprecision(7) << val;
                encoded.push_back(stream.str());
            };
            counter++;
        }
        endLoop:;
    }
    return encoded;
}

static vector<string> PSSMAAC(const string& seq, const string& seqName, vector<string> keys, const string& orderString, ifstream &file)
{
    vector<string> encoded;
    encoded.push_back(seqName);

    map<char, double> count;
    int l = seq.length();
    for(char c : orderString)
    {
        count[c] = 0;
    }
    // Read PSSM file
    string line;
    getline(file, line); // Skip the first line
    getline(file, line); // Skip the second line
    getline(file, line); // Skip the third line
    while (getline(file, line))
    {
        int counter = 0;
        istringstream iss(line);
        for(string line; iss >> line && counter < 22;) 
        {
            if(counter == 1)
            {
                bool exists = false;
                for (char c : orderString)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
            }
            if (counter >= 2)
            {
                int val = stoi(line);
                count[orderString[counter - 2]] += val;
            };
            counter++;
        }
        endLoop:;
    }
    for(char c : orderString)
    {
        count[c] /= (l * 1.0);
        stringstream stream;
        stream << fixed << setprecision(7) << count[c];
        encoded.push_back(stream.str());
    }
    return encoded;
}

static vector<string> BiPSSM(const string& seq, const string& seqName, vector<string> keys, const string& orderString, int n, ifstream &file)
{
    vector<string> encoded;
    encoded.push_back(seqName);

    map<string, double> count;
    int l = seq.length();
    for(int i = 0; i < keys.size(); i++)
    {
        count[keys[i]] = 0;
    }
    // Read PSSM file
    string line;
    int lines = 0;
    queue<array<int, 20>> prevLines;
    getline(file, line); // Skip the first line
    getline(file, line); // Skip the second line
    getline(file, line); // Skip the third line
    while (getline(file, line))
    {
        array<int, 20> newLine;
        int counter = 0;
        istringstream iss(line);
        for(string line; iss >> line && counter < 22;) 
        {
            if(counter == 1)
            {
                bool exists = false;
                for (char c : orderString)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
            }
            if (counter >= 2)
            {
                int val = stoi(line);
                newLine[counter - 2] = val;
                if (lines >= n)
                {
                    char c1 = orderString[counter - 2];
                    for (int i = 0; i < 20; i++)
                    {
                        char c2 = orderString[i];
                        int val2 = prevLines.front()[i];
                        string key(1, c1);
                        key.push_back(c2);
                        count[key] += (val * val2);
                    }
                }
            };
            counter++;
            if (counter == 22)
            {
                prevLines.push(newLine);
                if(lines >= n)
                    prevLines.pop();
                lines++;
            }
        }
        endLoop:;
    }

    for (string key : keys) 
    {
        count[key] /= (l - (n * 1.0));
        stringstream stream;
        stream << fixed << setprecision(7) << count[key];
        encoded.push_back(stream.str());
    }
    return encoded;
}

static vector<string> PSSMAC(const string& seq, const string& seqName, vector<string> keys, const string& orderString, int n, ifstream &file)
{
    vector<string> encoded;
    encoded.push_back(seqName);

    map<char, double> avg;
    map<string, double> count;
    int l = seq.length();
    for(char c : orderString)
    {
        avg[c] = 0;
    }
    for(string key : keys)
    {
        count[key] = 0;
    }
    // Read PSSM file
    string line;
    vector<array<int, 20>> lines;
    getline(file, line); // Skip the first line
    getline(file, line); // Skip the second line
    getline(file, line); // Skip the third line
    while (getline(file, line))
    {
        array<int, 20> newLine;
        int counter = 0;
        istringstream iss(line);
        for(string line; iss >> line && counter < 22;) 
        {
            if(counter == 1)
            {
                bool exists = false;
                for (char c : orderString)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
            }
            if (counter >= 2)
            {
                int val = stoi(line);
                newLine[counter - 2] = val;
                avg[orderString[counter - 2]] += val;
            };
            counter++;
            if (counter == 22)
            {
                lines.push_back(newLine);
            }
        }
        endLoop:;
    }

    for (char c : orderString)
    {
        avg[c] /= l;
    }

    for (int i = 1; i <= n; i++)
    {
        for (int j = 0; j < lines.size() - i; j++)
        {
            for (int k = 0; k < 20; k++)
            {
                char c = orderString[k];
                string key = to_string(i);
                key.push_back(c);
                count[key] += ((lines[j][k] - avg[c]) * (lines[j + i][k] - avg[c]));
            }
        }
    }
    
    for (string key : keys) 
    {
        count[key] /= (l - (n * 1.0));
        stringstream stream;
        stream << fixed << setprecision(7) << count[key];
        encoded.push_back(stream.str());
    }
    return encoded;
}

static vector<string> PPSSM(const string& seq, const string& seqName, vector<string> keys, const string& orderString, int n, ifstream &file)
{
    vector<string> encoded;
    encoded.push_back(seqName);
    vector<double> avg;
    map<char, double> avgChar;
    map<string, double> count;
    int l = seq.length();
    for(int i = 0; i < seq.size(); i++)
    {
        avg.push_back(0);
    }
    for(char c : orderString)
    {
        avgChar[c] = 0;
    }
    for(string key : keys)
    {
        count[key] = 0;
    }
    // Read PSSM file
    string line;
    vector<array<double, 20>> lines;
    getline(file, line); // Skip the first line
    getline(file, line); // Skip the second line
    getline(file, line); // Skip the third line
    int lineCount = 0;
    while (getline(file, line))
    {
        array<double, 20> newLine;
        int counter = 0;
        istringstream iss(line);
        for(string line; iss >> line && counter < 22;) 
        {
            if(counter == 1)
            {
                bool exists = false;
                for (char c : orderString)
                {
                    if (line[0] == c)
                    {
                        exists = true;
                        break;
                    }
                }
                if (!exists)
                    goto endLoop;
            }
            if (counter >= 2)
            {
                int val = stoi(line);
                newLine[counter - 2] = val;
                avg[lineCount] += val;
            };
            counter++;
            if (counter == 22)
            {
                lineCount++;
                lines.push_back(newLine);
            }
        }
        endLoop:;
    }

    for (int i = 0; i < l; i++)
    {
        avg[i] /= 20;
        double den = 0;
        for (int j = 0; j < 20; j++)
        {
            den += pow(lines[i][j] - avg[i], 2);
        }
        den = sqrt(den / 20);
        for (int j = 0; j < 20; j++)
        {
            lines[i][j] -= avg[i];
            lines[i][j] /= den;
            string key(1, orderString[j]);
            count["M" + key] += lines[i][j];
            if (i >= n)
                count["G" + key] += pow(lines[i - n][j] - lines[i][j], 2);
        }
    }

    int counter = 0;
    for (string key : keys) 
    {
        stringstream stream;
        if (counter < 20)
            count[key] /= l;
        else
            count[key] /= (l - n);
        stream << fixed << setprecision(7) << count[key];
        encoded.push_back(stream.str());
        counter++;
    }
    
    return encoded;
}

static vector<string> PseKRAAC(const string& seq, const string& seqName, const string& allowed, vector<string> keys, map<char, int> aaMap, const string& ktuple, 
    const string& subtype, int gapLambda)
{
    vector<string> encoded;
    encoded.push_back(seqName);

    int increment = subtype.compare("g-gap") == 0 ? gapLambda + 1 : 1;
    int rangeCheck = subtype.compare("g-gap") == 0 ? 1 : gapLambda;

    map<string, double> counts;

    if (ktuple == "1")
    {
        for (int i = 0; i < seq.length(); i += increment)
        {
            counts[to_string(aaMap[seq[i]] + 1)]++;
        }
    }
    else if (ktuple == "2")
    {
        for (int i = 0; i < seq.length(); i += increment)
        {
            if (i + rangeCheck < seq.length())
            {
                counts[to_string(aaMap[seq[i]] + 1) + "_" + to_string(aaMap[seq[i + rangeCheck]] + 1)]++;
            }
        }
    }
    else if (ktuple == "3")
    {
        for (int i = 0; i < seq.length(); i += increment)
        {
            if (i + rangeCheck < seq.length() && i + (2 * rangeCheck) < seq.length())
            {
                counts[to_string(aaMap[seq[i]] + 1) + "_" + to_string(aaMap[seq[i + rangeCheck]] + 1) + 
                    to_string(aaMap[seq[i + (2 * rangeCheck)]] + 1)]++;
            }
        }
    }

    for(int i = 0; i < keys.size(); i++) 
    {
        stringstream stream;
        stream << fixed << setprecision(7) << counts[keys[i]];
        encoded.push_back(stream.str());
    }
    return encoded;
}