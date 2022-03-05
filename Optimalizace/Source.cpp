//Simplex solver Spendley-Hext-Himsworth method
//Autor: Filip Nejedlý
//Git: https://github.com/Jacckii

//Libraries used:
//ImGui: 
//  Autor: Ocornut
//  Git: https://github.com/ocornut/imgui
//  Description: GUI Framework
//ImPlot:
//  Autor: epezent
//  Git: https://github.com/epezent/implo
//  Description: Plot widgets for ImGui

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <string> 
#include <vector>
#include <iostream>
#include <gui.h>
#include <implot.h>
#include <chrono>
#include "structs.h"

Point start(0.f, 0.f);
float delta = 0.2f;
float min_delta = 0.001f;
int max_iter_num = 300;
#define n 2

//Input function
float function(float x, float y) {
    float out = powf(1.5 - x + x * y, 2) + powf(2.25 - x + (x) * powf(y, 2), 2) + powf(2.625 - x + x * powf(y, 3), 2); //Beale function
    return out;
}

float function(Point a) {
    float x = a.x;
    float y = a.y;

    return function(x, y);
}

//Euclids ones vector
std::vector<std::vector<float>> euclid_vec;

//records of every calculated simplex value
std::vector<SimplexRecord> simplex_records;

/**
 * Generates simplex vector array
 * size - size of simplex vector array 
 */ 
void generateEuclidOneVec(int size) {
    for (int i = 0; i < size; i++) {
        std::vector<float> vec;

        for (int j = 0; j < size; j++){
            vec.push_back(j == i ? 1 : 0);
        }
        euclid_vec.push_back(vec);
    }
}

void printPoint(const char* name, Point a) {
    printf("Point %s : x: %.4f y: %.4f F(x): %.4f\n", name, a.x, a.y, function(a));
}

void printSimplex2D(Simplex2D& simplex) {
    printPoint("X1", simplex.x1);
    printPoint("X2", simplex.x2);
    printPoint("X3", simplex.x3);
}

/**
  * Generates starting points for simplex
  * ret - object containing generated simplex points
  */
Simplex2D generateStartingPoints() { 
    float d = delta * ((sqrt((float)n + 1.f) - 1.f) / ((float)n * sqrt(2.f)));
    float l = delta / sqrtf(2.f);

    Simplex2D simplex;
    simplex.x1 = start;

    simplex.x2.x = start.x + d * 1.f + l * euclid_vec[0][0];
    simplex.x2.y = start.y + d * 1.f + l * euclid_vec[0][1];

    simplex.x3.x = start.x + d * 1.f + l * euclid_vec[1][0];
    simplex.x3.y = start.x + d * 1.f + l * euclid_vec[1][1];

    return simplex;
}

/**
  * Finds the worst point, or second or n worst point
  * ignore: should we ignore some value?
  *  0 - Don't ignore
  *  1 - Ignore x1
  *  2 - Ignore x2
  *  3 - Ignore x3
  * returns:
  *  1 for x1,
  *  2 for x2,
  *  3 for x3 being the worst point
  */ 
int GetTheWorstPoint(Simplex2D points, int ignore = 0){
    float x1, x2, x3;

    x1 = function(points.x1);
    x2 = function(points.x2);
    x3 = function(points.x3);

    if (!ignore) {
        if (x1 > x2) {
            return x1 > x3 ? 1 : 3;
        }
        else {
            return x2 > x3 ? 2 : 3;
        }
    }
    else if (ignore == 1) {
        return x2 > x3 ? 2 : 3;
    }
    else if (ignore == 2) {
        return x1 > x3 ? 1 : 3;
    }
    else if (ignore == 3) {
        return x1 > x2 ? 1 : 2;
    }

    return 0;
}

/**
 * Calculates center between 2 points Xc
 * ret - returns Xc
 */
Point calculateCenter(Point a, Point b) {
    return (a + b) * (1.f / (float)n);
}

/**
  * point for pivoting:
  *  1 - x1
  *  2 - x2
  *  3 - x3
  * callculates Xp
  */
Simplex2D pivot(Simplex2D points, int point) {

    switch (point)
    {
    case 1:
        points.x1 = calculateCenter(points.x2, points.x3) * 2.f - points.x1; 
        break;
    case 2:
        points.x2 = calculateCenter(points.x1, points.x3) * 2.f - points.x2;
        break;
    case 3:
        points.x3 = calculateCenter(points.x1, points.x2) * 2.f - points.x3;
        break;
    default:
        std::runtime_error(std::string("Invalid parameter point:") + std::to_string(point) + " in pivot func!");
        break;
    }
    return points;
}

/**
 * Reduces simplex
 * ret - reduced simplex value
 */
Simplex2D reduce(Simplex2D& points, int Xq_index) {
    Simplex2D out;
    Point Xq;

    switch (Xq_index)
    {
    case 0:
        Xq = points.x1;
        break;
    case 1:
        Xq = points.x2;
        break;
    case 2:
        Xq = points.x3;
        break;
    default:
        std::runtime_error(std::string("Invalid parameter Xq_index:") + std::to_string(Xq_index) + " in reduce func!");
        break;
    }

    delta = delta / 2;
    out = points;
    
    if (Xq_index != 0)
        out.x1 = Xq + ((points.x1 - Xq) * 0.5f);

    if (Xq_index != 1)
        out.x2 = Xq + ((points.x2 - Xq) * 0.5f);

    if (Xq_index != 2)
        out.x3 = Xq + ((points.x3 - Xq) * 0.5f);

    return out;
}

/***
 * Simplex method (Spendley-Hext-Himsworth)
 */
void CalculateSimplex() {
    SimplexRecord rec;
    std::vector<int> last_iter;
    last_iter.assign(n + 1, 0);

    simplex_records.clear();
    generateEuclidOneVec(n);
    auto points = generateStartingPoints();
    int worst = 0, pre_worts = 0;
    int time_to_reduce = 4 * n;

    printSimplex2D(points);

    rec.simplex = points;
    rec.delta = delta;
    rec.iteration_index = 0;

    simplex_records.push_back(rec);
    for (int i = 1;; i++) {
        rec.simplex = points;
        rec.delta = delta;
        rec.iteration_index = i;

        simplex_records.push_back(rec);

        if (delta <= min_delta || i >= max_iter_num) {
            break;
        }

        bool should_continue = false;
        for (unsigned int j = 0; j < last_iter.size(); j++) {
            if (last_iter[j] + time_to_reduce > i)
                continue;

            points = reduce(points, j);
            std::cout << "Simplex reduced around point x" << std::to_string(j + 1) << std::endl;

            printSimplex2D(points);

            for (int g = 0; g < n + 1; g++) {
                last_iter[g] = i;
            }

            worst = 0;
            should_continue = true;
        }

        if (should_continue)
            continue;

        pre_worts = worst;
        worst = GetTheWorstPoint(points);

        if (worst == pre_worts) {
            worst = GetTheWorstPoint(points, pre_worts);
        }

        std::cout << "Point for pivoting is x" << worst << std::endl;

        points = pivot(points, worst);
        last_iter[worst - 1] = i;

        printSimplex2D(points);
    }

    std::cout << "Result is:" << std::endl;
    printSimplex2D(points);
}

// ------------------------- MAIN, GUI -----------------------------------------

int main(int argc, char* argv[]) {
	std::unique_ptr<GUI> gui = std::make_unique<GUI>(L"2D Simplex solver", 900, 580);
	try {
		while (true) {
			if (gui->beginFrame()) {
                static ImGuiWindowFlags flags = ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoCollapse;
                const ImGuiViewport* viewport = ImGui::GetMainViewport();
                ImGui::SetNextWindowPos(viewport->WorkPos);
                ImGui::SetNextWindowSize(viewport->WorkSize);
				ImGui::Begin("Simplex solver (Spendley-Hext-Himsworth)", NULL, flags);

                ImGui::Columns(2, NULL, false);
                ImGui::SetColumnWidth(0, 580);
                {
                    static bool show_min = true;
                    static int pos = 0;
                    static float mk_size = ImPlot::GetStyle().MarkerSize;
                    static float mk_weight = ImPlot::GetStyle().MarkerWeight;
                    static float scale_min = 999999.f;
                    static float scale_max = 200.f;
                    static Point min_pos;
                    ImPlot::PushColormap(ImPlotColormap_Viridis);
                    ImPlot::ColormapScale("##HeatScale", scale_min, scale_max, ImVec2(60, 500));
                    ImGui::SameLine();

                    int resolution = 200;
                    static std::vector<float> data;
                    if (data.size() == 0) {
                        data.resize(resolution * resolution);
                        float step = 8.f / (float)resolution;
                        for (int i = 0; i < resolution; i++) {
                            float y = 4.f - i * step;
                            for (int j = 0; j < resolution; j++) {
                                float x = -4.f + j * step;
                                auto out = function(x, y);
                                data[i * resolution + j] = out;

                                /*if (out > scale_max)
                                    scale_max = out;*/

                                if (out < scale_min) {
                                    scale_min = out;
                                    min_pos.x = x;
                                    min_pos.y = y;
                                }
                            }
                        }
                    }

                    if (ImPlot::BeginPlot("##Heatmap1", ImVec2(500, 500))) {
                        ImPlot::SetupAxes("X", "Y", NULL, NULL);
                        ImPlot::SetupAxesLimits(-4, 4, -4, 4);
                        ImPlot::PlotHeatmap("##T", data.data(), resolution, resolution, scale_min, scale_max, NULL, ImPlotPoint(-4, -4), ImPlotPoint(4, 4));
                        
                        if (show_min) {
                            auto circle_size = 2.5f;
                            ImVec2 cntr = ImPlot::PlotToPixels(ImPlotPoint(min_pos.x, min_pos.y));
                            ImPlot::GetPlotDrawList()->AddCircleFilled(cntr, circle_size, IM_COL32(200, 50, 20, 255), 20);
                        }

                        if (simplex_records.size() > pos) {
                            auto simplex = simplex_records[pos].simplex;
                            ImVec2 pointx1 = ImPlot::PlotToPixels(ImPlotPoint(simplex.x1.x, simplex.x1.y));
                            ImVec2 pointx2 = ImPlot::PlotToPixels(ImPlotPoint(simplex.x2.x, simplex.x2.y));
                            ImVec2 pointx3 = ImPlot::PlotToPixels(ImPlotPoint(simplex.x3.x, simplex.x3.y));
                            ImPlot::GetPlotDrawList()->AddLine(pointx1, pointx2, IM_COL32(200, 50, 20, 255), 1.f);
                            ImPlot::GetPlotDrawList()->AddLine(pointx2, pointx3, IM_COL32(200, 50, 20, 255), 1.f);
                            ImPlot::GetPlotDrawList()->AddLine(pointx3, pointx1, IM_COL32(200, 50, 20, 255), 1.f);
                        }
                        ImPlot::EndPlot();
                    }
                    ImPlot::PopColormap();

                    ImGui::NextColumn(); 
                    {
                        static float start_delta = delta;
                        ImGui::Text("Plot settings:");
                        ImGui::SliderFloat("Min", &scale_min, 0.f, 400.f);
                        ImGui::SliderFloat("Max", &scale_max, 0.f, 400.f);
                        ImGui::Checkbox("Show minimum point", &show_min);
                        ImGui::Text("Simplex settings:");
                        ImGui::SliderFloat("Start x", &start.x, -5, 5);
                        ImGui::SliderFloat("Start y", &start.y, -5, 5);
                        ImGui::SliderFloat("Start delta", &start_delta, 0.f, 1.f);
                        ImGui::SliderFloat("Min delta", &min_delta, 0.f, 0.1f, "%.6f");
                        ImGui::SliderInt("Max iterace", &max_iter_num, 0, 1000);
                        ImGui::Text("Calculations:");
                        if (ImGui::Button("Calculate!", ImVec2(-1.f, 0.f))) {
                            delta = start_delta;
                            CalculateSimplex();
                            pos = 0;
                        }
                        if (simplex_records.size() > 0) {
                            static bool playing = false;
                            static std::chrono::system_clock::time_point last_time_stamp = std::chrono::system_clock::now();
                            auto step = 80;
                            auto max_ = simplex_records.size() - 1;
                            ImGui::Text("%d/%d Krok iterace", pos, max_);
                            ImGui::SliderInt("##timeline", &pos, 0, max_);
                            ImGui::SameLine();
                            if (ImGui::Button(playing ? "Stop!" : "Play!", ImVec2(-1.f, 0.f))) {
                                if (!playing) {
                                    last_time_stamp = std::chrono::system_clock::now();
                                    pos = 0;
                                    playing = true;
                                }
                                else {
                                    playing = false;
                                }
                            }

                            if (playing) {
                                double elapsed_time_ms = std::chrono::duration<double, std::milli>(std::chrono::system_clock::now() - last_time_stamp).count();
                                if (elapsed_time_ms > step) {
                                    pos++;
                                    last_time_stamp = std::chrono::system_clock::now();
                                    if (pos >= max_)
                                        playing = false;
                                }
                            }

                            ImGui::NewLine();
                            if (simplex_records.size() > 0 && simplex_records.size() > pos) {
                                auto record = simplex_records[pos];
                                ImGui::Text("Simplex iteration info:\n - iteration number: %d\n - Pos:\n - delta: %.6f"
                                    "\n\t\tX1: (%.2f, %.2f) f(X1) = %.2f\n\t\tX2: (%.2f, %.2f) f(X2) = %.2f\n\t\tX3: (%.2f, %.2f) f(X3) = %.2f", record.iteration_index, record.delta,
                                    record.simplex.x1.x, record.simplex.x1.y, function(record.simplex.x1), record.simplex.x2.x, record.simplex.x2.y, 
                                    function(record.simplex.x2),record.simplex.x3.x, record.simplex.x3.y, function(record.simplex.x3));

                                auto last = simplex_records[simplex_records.size() - 1];
                                ImGui::Text("Solution is: %.2f, %.2f", last.simplex.x1.x,
                                    last.simplex.x1.y);
                            }
                            else {
                                ImGui::Text("Simplex iteration info:\n");
                            }
                        }
                        ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.f, 1.f, 1.f, 0.5f));
                        ImGui::TextWrapped("To change calculated function please edit Source.cpp line 36!");
                        ImGui::PopStyleColor();
                    }
                }
				ImGui::End();
				gui->endFrame();
			}

			if (gui->getMsg().message == WM_QUIT) {
				break;
			}
		}
	}
	catch (const std::exception& ex) {
        printf("\n\n!!!\n%s\n!!\n", ex.what()); 
	}
    return 0;
}