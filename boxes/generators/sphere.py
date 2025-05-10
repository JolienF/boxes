# Copyright (C) 2013-2016 Jolien Franke
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

# A gore is term in globe-making that describes a lens-shaped figure (with a map usually) that is wrapped around a sphere to create a globe
# The coordinate system of the 3D object is x,y,z. As to not confuse it with the coordinate system of the gores, we use x,u (x is the same in both systems)
# Explanation of the math will be in a separate file

import math

import numpy
from boxes import *
from collections import namedtuple


class Sphere(Boxes):
    """Actually not a sphere, but a hosohedron. Also not actually a box, but a globe, lamp or ornament"""


# TODO thickness/ 2 must be larger than kerf or other way around
# Tabs to be even number! (and at least 2)
    ui_group = "Misc"

    def __init__(self) -> None:
        super().__init__()
        self.argparser.add_argument(
            "--sphere_radius",  action="store", type=float, default=200,
            help="The radius of the assembled sphere")
        self.argparser.add_argument(
            "--amount_gores",  action="store", type=int, default=6, #TODO minium = 3
            help="The amount of gores/parts you want the sphere to have, has to be at least 3")
        self.argparser.add_argument(
            "--top_hole_radius",  action="store", type=float, default=30,
            help="The size of the circular hole at the top") #TODO can't be negative, test for 0
        self.argparser.add_argument(
            "--bottom_hole_radius", action="store", type=float, default=120,
            help="The size of the polygonal hole at the bottom")  #TODO can't be negative, test for 0


        defaultgroup = self.argparser._action_groups[1]
        for action in defaultgroup._actions:
            if action.dest == 'tabs':
                action.type = int
                action.default = 10
                action.help = "number of tabs"
            if action.dest == 'thickness':
                action.default = 1.0
        defaultgroup.add_argument( # I placed it here, hoping it would be grouped together with tabs, but this doesn't work (tips welcome)
            "--corner_tab", action="store", type=float, default=4.0,
            help="The length of the tabs on the corners (in mm)(not supported everywhere). Keep as small as your material strength allows for cleaner result")
        self.argparser.add_argument(
            "--tab_width", action="store", type=float, default=8.0, # has to be larger than material thickness
            help="The width of the tabs (in mm)")
        # self.argparser.add_argument(
        #     "--thickness", action="store", type=float, default=1.0,
        #     help="The thickness of the material (in mm)")

        # parse = argparse.ArgumentParser(conflict_handler="resolve")
        # parse.add_argument(
        #     "--tabs", action="store", type=float, default=0.0,
        #     help="new tabs")

    Curve = namedtuple('Curve', ["degrees", "radius"])


    def calculateXOfGore(self, u):

        return math.sin(u / self.sphere_radius) * self.halfBellyLens

    def coordinatesPartOfGore(self, u_start, u_stop, isRight):
        N = self.resolution
        points = []

        direction = 1
        if (not isRight):
            direction = -1

        for i in range (N + 1):
            u = (u_stop - u_start) / N * i + u_start
            x = direction * self.calculateXOfGore(u)
            points.append((x, u))

        # self.drawPoints(points, 1, close=False) #kerfdir = 1

        return(points)

    def coordinatesPartOfOffset(self, u_start, u_stop, normalDistance, isRight):                           # see wiki of Parallel Curves x_d(t) = x(t) + d * n(t)
        N = self.resolution
        points = []

        direction = 1
        if (not isRight):
            direction = -1

        for i in range (N + 1):
            u = (u_stop - u_start) / N * i + u_start
            u_offset = u + normalDistance * math.sin(self.calculateNormalAngle(u))
            x_offset = direction * (self.calculateXOfGore(u) + normalDistance * math.cos(self.calculateNormalAngle(u))) # could be replaced by self.calculateXofOffsetGore
            points.append((x_offset, u_offset))

        # self.drawPoints(points, 1, close=False)
        return(points)



    def calculateNormalAngle(self, u):
        return self.calculateTangentAngle(u) - 90 * math.pi / 180

    def calculateTangentAngle(self, u):                                                                              #derivatives of u (goreHeigth) and x (cos(pi * i) * a * pi)
        return math.atan2(self.gore_heigth, math.cos(math.pi * (u / self.gore_heigth)) * math.pi * self.halfBellyLens) #atan2 to prevent division by 0 and quadrant (opposite sides not possible now)

    def mapYToU(self, y):
        return (y / math.pi) * self.gore_heigth

    def normalCompensation(self, u):                                            #So the offset can be drawn up to the same horizonal line as the equivalent u of the gore
        return math.sin(self.calculateNormalAngle(u)) * self.tab_width

    def calculateXofTopAndGoreIntersection(self):
        return math.sqrt(self.halfBellyLens**2 * self.top_hole_radius**2 / (self.halfBellyLens**2 + self.sphere_radius**2))

    def calculateXofOffsetGore(self, u):
        return self.calculateXOfGore(u) + math.cos(self.calculateNormalAngle(u)) * self.tab_width

    def calculateUpperUOfGore(self, x):
        return (math.pi - math.asin(x / self.halfBellyLens)) * self.sphere_radius

    def calculateUOfBottomHole(self):
        y = math.asin(self.bottom_hole_radius / self.sphere_radius)
        return self.mapYToU(y)

    def plotBottomHoleByRadius(self):
        self.ctx.move_to(-self.x_rightGoreBottom, self.u_rightGoreBottom)
        self.ctx.line_to(self.x_rightGoreBottom, self.u_rightGoreBottom)

        # def plotBottomHoleByRadius(self):
        # self.medium.penup()
        # self.medium.setposition(-self.x_rightGoreBottom, self.u_rightGoreBottom)
        # self.medium.pendown()
        # self.medium.goto(self.x_rightGoreBottom, self.u_rightGoreBottom)

    def coordinatesTopHole(self, x_start, x_stop):
        N = self.resolution
        points = []

        for i in range (N + 1):
            x = (x_stop - x_start) / N * i + x_start
            y = math.acos(-(math.sqrt(self.sphere_radius**2 - self.top_hole_radius**2 + x**2)) / self.sphere_radius)
            u = self.mapYToU(y)
            points.append((x, u))

        return points

    def coordinatesToPolyline(self, points):
        # needs all the coordinates at once, if this could handle separate pieces, I would need to remember the absolute angle
        polyPoints = [0]
        distance = 0
        previousWasCurve = False
        previousAngle = 0
        i = 1

        while i < len(points):
            if (isinstance(points[i], self.Curve)):
                polyPoints.extend([points[i], points[i + 1]])
                previousAngle += points[i].degrees
                previousWasCurve = True
                i += 2

            else:
                if (previousWasCurve):
                    previousWasCurve = False
                    i += 1
                    continue

                absoluteAngle = math.degrees(math.atan2(points[i][1] - points[i - 1][1], points[i][0] - points[i - 1][0]))
                relativeAngle = absoluteAngle - previousAngle
                relativeAngle = (relativeAngle + 180) % 360 - 180

                distance = math.dist(points[i], points[i - 1])

                if (distance > 0.001):
                    polyPoints.append(relativeAngle)
                    polyPoints.append(distance)
                    previousWasCurve = False
                    previousAngle = absoluteAngle
                i += 1

        return(polyPoints)



        # Tabs
        # Surrounding spaces don't make much sense in this case, so I've left those out
        # I've changed tabs from mm to amount, because the first and last tab are very important to keep the shape nice

    def coordinatesTabStart(self, isRight):
        direction = 1
        if (not isRight):
            direction = -1

        return [self.Curve(-180 * direction, self.thickness / 2), 0, self.Curve(180 * direction, (self.tab_width - self.thickness) / 2), 0]

    def coordinatesTabEnd(self, isRight):
        direction = 1
        if (not isRight):
            direction = -1

        return [self.Curve(180 * direction, (self.tab_width - self.thickness) / 2), 0, self.Curve(-180 * direction, self.thickness / 2), 0]

    # @holeCol
    # def plotTabStart(self, u, isRight):
    #     direction = 1
    #     if (not isRight):
    #         direction = -1

    #     # self.ctx.line_to(self.calculateXOfGore(u), u)
    #     # self.moveTo(startPoints[0], startPoints[1])
    #     self.moveTo(self.calculateXOfGore(u) * direction, u) #TODO firstfinger
    #     # self.ctx.move_to(self.calculateXOfGore(u), u)
    #     self.ctx.rotate(direction * (self.calculateTangentAngle(u)))

    #     # self.left(direction * (self.calculateTangentAngle(u) * 360 / (2* math.pi) - 180))
    #     self.corner(-180, self.thickness / 2)

    #     self.corner(180, (self.tab_width - self.thickness) / 2)


        # self.medium.setheading(0)
        # self.medium.left(direction * (self.calculateTangentAngle(u) * 360 / (2* math.pi) - 180))
        # self.medium.circle(self.materialRadius, -180)
        # self.medium.right(180)
        # self.medium.circle(self.fingerRadius, 180)

    # @holeCol
    # def plotTabEnd(self, u, isRight, startPoints):
    #     direction = 1
    #     if (not isRight):
    #         direction = -1

    #     self.moveTo(startPoints[0], startPoints[1]) #TODO firstfinger
    #     self.ctx.rotate(direction * (self.calculateTangentAngle(u)))


    #     self.corner(180, (self.tab_width - self.thickness) / 2)
    #     self.corner(-180, self.thickness / 2)


        # This function divides the tabs over the length of u, not the length of the gore itself. Shortcuts were taken
    def divideGore(self, numberOfTabs, tinyTabLength):
        return numpy.linspace(self.u_rightGoreBottom + tinyTabLength,
                              self.u_rightGoreTop - tinyTabLength,
                              num=numberOfTabs + 1)

    def render(self):
        tabs = self.tabs
        allCoordinates = []

        self.resolution = int(self.sphere_radius / 10)  # This is arbitrary. I just want the resolution to be proportional to the total size
        self.gore_heigth = math.pi * self.sphere_radius  # the midline of the gore is the same as the line on the sphere, so half a circumference
        self.bellyLens = 2 * math.tan((2 * math.pi) / (2 * self.amount_gores)) * self.sphere_radius
        self.halfBellyLens = self.bellyLens / 2

        self.x_rightGoreTop = self.calculateXofTopAndGoreIntersection() #TODO remove right, it's confusing
        self.u_rightGoreTop = self.calculateUpperUOfGore(self.x_rightGoreTop)
        self.u_rightGoreBottom = self.calculateUOfBottomHole()
        self.x_rightGoreBottom = self.calculateXOfGore(self.u_rightGoreBottom)

        self.u_tabPoints = self.divideGore(self.tabs, self.corner_tab)

        #tabcheck
        # for i in range(len(self.u_tabPoints)):
        #     self.ctx.move_to(0, self.u_tabPoints[i])
        #     self.ctx.line_to(150, self.u_tabPoints[i])
        # self.drawPoints(self.coordinatesPartOfGore(self.u_rightGoreBottom, self.u_rightGoreTop, True))

        allCoordinates.extend(self.coordinatesPartOfGore(self.u_rightGoreBottom, self.u_tabPoints[0], True))

        for i in range(0, len(self.u_tabPoints) - 1, 2):
            allCoordinates.extend(self.coordinatesTabStart(True))
            allCoordinates.extend(self.coordinatesPartOfOffset(self.u_tabPoints[i], self.u_tabPoints[i + 1], self.tab_width, True))
            allCoordinates.extend(self.coordinatesTabEnd(True))
            allCoordinates.extend(self.coordinatesPartOfGore(self.u_tabPoints[i + 1], self.u_tabPoints[i + 2], True))

        lastTabPoint = self.u_tabPoints[-1]
        upperCornerU = self.u_rightGoreTop - self.normalCompensation(self.u_tabPoints[-1])

        allCoordinates.extend(self.coordinatesTabStart(True))

        if (lastTabPoint < upperCornerU):
            allCoordinates.extend(self.coordinatesPartOfOffset(lastTabPoint, upperCornerU, self.tab_width, True)) # line is not entirely straight, because te normalcompensation is simplified (should be at a bit lower u)
        else:
            allCoordinates.extend([(self.calculateXofOffsetGore(upperCornerU), upperCornerU + self.tab_width * math.sin(self.calculateNormalAngle(upperCornerU)))])
        allCoordinates.extend([(self.x_rightGoreTop, self.u_rightGoreTop)])

        allCoordinates.extend(self.coordinatesTopHole(self.x_rightGoreTop, -self.x_rightGoreTop))

        #Left half
        allCoordinates.extend(self.coordinatesPartOfGore(self.u_rightGoreTop, lastTabPoint, False))

        for i in range(len(self.u_tabPoints) - 1, 0, -2):
            allCoordinates.extend(self.coordinatesTabStart(True))
            allCoordinates.extend(self.coordinatesPartOfOffset(self.u_tabPoints[i], self.u_tabPoints[i - 1], self.tab_width, False))
            allCoordinates.extend(self.coordinatesTabEnd(True))
            allCoordinates.extend(self.coordinatesPartOfGore(self.u_tabPoints[i - 1], self.u_tabPoints[i - 2], False))


        firstTabPoint = self.u_tabPoints[0]
        lowerCornerU = self.u_rightGoreBottom + self.normalCompensation(self.u_tabPoints[-1])
        # bottomTabPoints = []

        allCoordinates.extend(self.coordinatesTabStart(True))

        if (firstTabPoint > lowerCornerU):
            allCoordinates.extend(self.coordinatesPartOfOffset(firstTabPoint, lowerCornerU, self.tab_width, False))
        else:
            allCoordinates.extend([(-self.calculateXofOffsetGore(lowerCornerU), lowerCornerU + self.tab_width * math.sin(self.calculateNormalAngle(lowerCornerU)))])
        allCoordinates.extend([(-self.x_rightGoreBottom, self.u_rightGoreBottom)])

        allCoordinates.extend([(self.x_rightGoreBottom, self.u_rightGoreBottom)])

        polyPoints = self.coordinatesToPolyline(allCoordinates)



        self.moveTo(-self.halfBellyLens, 30)
        moveX = self.bellyLens + self.tab_width * 2 + 30
        for i in range(self.amount_gores):
            self.moveTo(moveX, 0)
            self.polyline(*polyPoints)




















        # thisPoint = self.coordinatesPartOfGore(self.u_rightGoreBottom, self.u_tabPoints[0], True)
        # # self.ctx.line_to(self.calculateXOfGore(self.u_tabPoints[1]), self.u_tabPoints[1])

        # for i in range(0, len(self.u_tabPoints) - 1, 2):
        #     # self.plotTabStart(self.u_tabPoints[i], True, color=Color.ANNOTATIONS)
        #     self.plotTabStart(self.u_tabPoints[i], True)
        #     punt = self.coordinatesPartOfOffset(self.u_tabPoints[i], self.u_tabPoints[i + 1], self.tab_width, True)
        #     self.plotTabEnd(self.u_tabPoints[i + 1], True, punt)
        #     thisPoint = self.coordinatesPartOfGore(self.u_tabPoints[i + 1], self.u_tabPoints[i + 2], True)

        # lastTabPoint = self.u_tabPoints[-1]
        # upperCornerU = self.u_rightGoreTop - self.normalCompensation(self.u_tabPoints[-1])

        # self.plotTabStart(lastTabPoint, True)

        # if (lastTabPoint < upperCornerU):
        #     self.coordinatesPartOfOffset(lastTabPoint, upperCornerU , self.tab_width, True)
        #     self.ctx.move_to(self.calculateXofOffsetGore(upperCornerU), self.u_rightGoreTop)
        # else:
        #     self.ctx.move_to(self.calculateXofOffsetGore(lastTabPoint), lastTabPoint + self.normalCompensation(lastTabPoint))
        # self.ctx.line_to(self.x_rightGoreTop, self.u_rightGoreTop)

        # self.coordinatesTopHole(self.x_rightGoreTop, -self.x_rightGoreTop)

        # self.coordinatesPartOfGore(self.u_rightGoreTop, self.u_tabPoints[-1], False)

        # for i in range(len(self.u_tabPoints) - 1, 0, -2):
        #     self.plotTabStart(self.u_tabPoints[i], False)
        #     punt = self.coordinatesPartOfOffset(self.u_tabPoints[i], self.u_tabPoints[i - 1], self.tab_width, False)
        #     self.plotTabEnd(self.u_tabPoints[i - 1], False, punt)
        #     thisPoint = self.coordinatesPartOfGore(self.u_tabPoints[i - 1], self.u_tabPoints[i - 2], False)

        # firstTabPoint = self.u_tabPoints[0]
        # lowerCornerU = self.u_rightGoreBottom + self.normalCompensation(self.u_tabPoints[-1])
        # bottomTabPoints = []
        # self.plotTabStart(firstTabPoint, False)

        # if (firstTabPoint > lowerCornerU):
        #     self.coordinatesPartOfOffset(self.u_tabPoints[0], lowerCornerU , self.tab_width, False)
        #     self.ctx.move_to(-self.calculateXofOffsetGore(lowerCornerU), self.u_rightGoreBottom)  # proper fix would be to make it one object, so I dont need to move to a point, would save some functions and calculations
        #     # bottomTabPoints.append((-self.calculateXofOffsetGore(lowerCornerU), self.u_rightGoreBottom))
        # else:
        #     self.ctx.move_to(-self.calculateXofOffsetGore(firstTabPoint), firstTabPoint + self.normalCompensation(firstTabPoint))
        #     # bottomTabPoints.append((-self.calculateXofOffsetGore(firstTabPoint), float(firstTabPoint + self.normalCompensation(firstTabPoint))))
        # self.ctx.line_to(-self.x_rightGoreBottom, self.u_rightGoreBottom)




        # bottomTabPoints.append((-self.x_rightGoreBottom, self.u_rightGoreBottom))
        # print(firstTabPoint)
        # print(self.normalCompensation(firstTabPoint))
        # print(bottomTabPoints)
        # print(firstTabPoint + self.normalCompensation(firstTabPoint))
        # self.drawPoints(bottomTabPoints, 1, False)
        # self.plotBottomHoleByRadius()

        # print(self.u_tabPoints)
        # self.edge(tabs)hnqvb
