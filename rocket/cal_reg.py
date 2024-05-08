'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2024-03-26 15:27:46
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2024-03-26 16:24:30
FilePath: /pulsar/rocket/cal_reg.py
Description: 
Copyright (c) 2024 by baotong, All Rights Reserved. 
'''
from shapely.geometry import Point
from shapely.ops import unary_union

def parse_region_file(file_path):
    """
    解析region文件，返回圆环和小圆的shapely几何对象列表。

    Args:
    file_path (str): region文件路径.

    Returns:
    tuple: 包含圆环和小圆shapely几何对象的元组.
    """
    annulus = None
    circles = []
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('(')
            shape_type = parts[0]
            coords = parts[1].replace(')', '').split(',')
            if shape_type == 'annulus':
                center_x, center_y, inner_radius, outer_radius = map(float, coords)
                annulus = Point(center_x, center_y).buffer(outer_radius).difference(Point(center_x, center_y).buffer(inner_radius))
            elif shape_type == '-circle':
                center_x, center_y, radius = map(float, coords)
                circles.append(Point(center_x, center_y).buffer(radius))
    # print(annulus,circles)
    return annulus, circles

def calculate_desired_area(annulus, circles):
    # 如果没有小圆，则返回整个圆环的面积
    if not circles:
        return annulus.area
    # 计算所有小圆的并集
    circles_union = unary_union(circles)
    # 计算小圆之间的相交部分
    circles_intersection = unary_union(circles).intersection(unary_union(circles))
    # 计算圆环与所有小圆的相交部分
    annulus_intersect_circles = annulus.intersection(circles_union)
    # 计算不在任何一个小圆内的面积
    desired_area = annulus.difference(circles_union).area
    return desired_area

path='/Volumes/pulsar/SgrA/merge_data/timing/region_startover/region_242/region_90/'
annulus, circles = parse_region_file(path+'1532_bkg.reg')
desired_area = calculate_desired_area(annulus, circles)

print("所需面积:", desired_area)




