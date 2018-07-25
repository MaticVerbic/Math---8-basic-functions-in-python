import math
import matplotlib.pyplot as plt

class Linear:

    def __init__(self, k=1, n=0):
        #k*x + n
        self.k  = k
        self.n = n


    def x_intercept(self):
        return -self.n/self.k

    def y_intercept(self):
        return self.n

    def angle(self):
        return self.k

    def coordinates(self, l=[]):
        return [(x, (self.k*x)+self.n) for x in l]

    def derivative(self):
        return n


    def draw(self, x):
        coords = self.coordinates(x)
        xs = [x for x, y in coords]
        ys = [y for x, y in coords]
        plt.ylim(-20, 20)
        plt.xlim(-20, 20)
        plt.plot(xs, ys, "-r")
        plt.show()

    def __add__(self, other):
        return Linear(self.k+other.k, self.n+other.n)

    def __sub__(self, other):
        return Linear(self.k-other.k, self.n-other.n)

    def __mul__(self, other):
        return Linear(self.k * other.k, self.n * other.n)

    def __truediv__(self, other):
        return Linear(self.k / other.k, self.n / other.n)

    def __repr__(self):
        return "f(x) = {}x{}{}".format(self.k , "+" if self.n >= 0 else "-", abs(self.n ))

    def __str__(self):
        return "f(x) = {}x{}{}".format(self.k, "+" if self.n > -1 else "-", abs(self.n))


class Quadratic:
    def __init__(self, a, b, c):
        #a*x**2 + b*x + c
        self.a = a
        self.b = b
        self.c = c

    def __repr__(self):
        return "f(x) = {}x^2{}{}x{}{}".format(self.a, "+" if self.b >= 0 else "-", self.b,
                                                "+" if self.c >= 0 else "-", self.c)

    def __str__(self):
        return "f(x) = {}x^2{}{}x{}{}".format(self.a, "+" if self.b >= 0 else "-", abs(self.b),
                                                "+" if self.c >= 0 else "-", abs(self.c))

    def __add__(self, other):
        return Quadratic(self.a+other.a, self.b+other.b, self.c+other.c)

    def __sub__(self, other):
        return Quadratic(self.a - other.a, self.b - other.b, self.c - other.c)

    def __mul__(self, other):
        return Quadratic(self.a * other.a, self.b * other.b, self.c * other.c)

    def __truediv__(self, other):
        return Quadratic(round(self.a / other.a, 2), round(self.b / other.b, 2), round(self.c / other.c, 2))

    def discriminant(self):
        return self.b**2 - (4*self.a*self.c)

    def x_intercept(self):
        return [(round((-self.b - math.sqrt(self.discriminant()))/(2*self.a), 2), 0),
                (round((-self.b + math.sqrt(self.discriminant()))/(2*self.a), 2), 0)]

    def y_intercept(self):
        return self.c

    def vertex(self):
        x = -self.b/(2*self.a)
        y = self.a*(x**2) + self.b*x + self.c
        return (x, y)

    def y(self, x):
        return (self.a*(x**2)) + (self.b * x) + self.c

    def coordinates(self, l):
        return [(x, self.y(x)) for x in l]

    def derivative(self):
        return Quadratic(0, 2*self.a, self.b)


    def draw(self, x):
        coords = self.coordinates(x)
        xs = [x for x, y in coords]
        ys = [y for x, y in coords]
        plt.ylim(-20, 20)
        plt.xlim(-20, 20)
        plt.plot(xs, ys, "-r")
        plt.show()




class Exponential:
    def __init__(self, a=1, b=1, c=1 , d = 0):
        #a * (b**(x+c)) + d
        self.a = a
        self.b = b
        self.c = c
        self.d = d

    def y_intercept(self):
        return (0, self.a * (self.b**(0+self.c)) + self.d)

    def horizontal_asymptote(self):
        if self.d != 0:
            return self.d
        else:
            return None

    def coordinates(self, l):
        return [(x, round(self.a * self.b**(x+self.c) + self.d, 2)) for x in l]

    def draw(self, x):
        coords = self.coordinates(x)
        xs = [x for x, y in coords]
        ys = [y for x, y in coords]
        plt.ylim(-20, 20)
        plt.xlim(-20, 20)
        plt.plot(xs, ys, "-r")
        plt.show()

    def __repr__(self):
        return "f(x): {}*({}^x)".format(self.a, self.b) if self.d == 0 \
            else "f(x): {}*({}^x){}{}".format(self.a, self.b, "+" if self.d>0 else "-", abs(self.d))

    def __str__(self):
        return "f(x): {}*({}^x)".format(self.a, self.b) if self.d == 0 \
            else "f(x): {}*({}^x){}{}".format(self.a, self.b, "+" if self.d>0 else "-", abs(self.d))


class Errors(Exception):
    pass


class Polynomial:
    def __init__(self, constants, degrees):
        #constants = c, degrees = d (c[0]*x^d[0] +- c[1]*x^d[1] +- ... +- c[n]*x^d[n],
        #degree of + constant = 0 or no-degree
        self.constants = constants
        self.degrees = degrees

        try:
            error = Errors
            if (len(self.degrees) not in range(len(self.constants)-1, len(self.constants)+1)):
                raise error
        except error:
            print("Polynomial error: Length of constants must be equal or 1 larger than length of degrees")
            exit()
        finally:
            if len(self.degrees) == len(self.constants):
                self.expression = [(const, deg) for const, deg in zip(self.constants, self.degrees)]
            else:
                self.expression = [(const, deg) for const, deg in zip(self.constants, self.degrees)] \
                                  + [(self.constants[-1], 0)]
            self.expression.sort(key=lambda x: x[1])
            self.expression = self.expression[::-1]

    def __repr__(self):
        return "{}x^{}".format(self.expression[0][0], self.expression[0][1]) + \
            "".join(["{}{}{}{}{}".format(
            "+" if constant >= 0 else "", constant, "x" if degree > 1 else "",  "^" if degree > 1 else "" ,degree if degree > 1 else "")
            for constant, degree in self.expression[1:]])

    def __str__(self):
        return "{}x^{}".format(self.expression[0][0], self.expression[0][1]) + \
            "".join(["{}{}{}{}{}".format(
            "+" if constant >= 0 else "", constant, "x" if degree > 1 else "",  "^" if degree > 1 else "" ,degree if degree > 1 else "")
            for constant, degree in self.expression[1:]])

    def __add__(self, other):
        constants = []
        degrees = []
        for constant, degree in self.expression:
            for constant_o, degree_o in other.expression:
                if degree == degree_o and degree not in degrees:
                    constants.append(constant + constant_o)
                    degrees.append(degree)
                    break
                else:
                    continue

        for constant, degree in self.expression:
            if degree not in degrees:
                constants.append(constant)
                degrees.append(degree)

        for constant, degree in other.expression:
            if degree not in degrees:
                constants.append(constant)
                degrees.append(degree)

        i = 0
        while i < len(constants)-1:
            if constants[i] == 0:
                del constants[i]
                del degrees[i]
                continue
            else:
                i += 1

        return Polynomial(constants, degrees)

    def __sub__(self, other):
        constants = []
        degrees = []
        for constant, degree in self.expression:
            for constant_o, degree_o in other.expression:
                if degree == degree_o and degree not in degrees:
                    constants.append(constant - constant_o)
                    degrees.append(degree)
                    break
                else:
                    continue

        for constant, degree in self.expression:
            if degree not in degrees:
                constants.append(constant)
                degrees.append(degree)

        for constant, degree in other.expression:
            if degree not in degrees:
                constants.append(constant)
                degrees.append(degree)

        i = 0
        while i < len(constants) - 1:
            if constants[i] == 0:
                del constants[i]
                del degrees[i]
                continue
            else:
                i += 1

        return Polynomial(constants, degrees)

    def __mul__(self, other):
        constants = []
        degrees = []
        muls = []
        for constant, degree in self.expression:
            for other_constant, other_degree in other.expression:
                muls.append((constant*other_constant, degree+other_degree))

        for ind, (constant, degree) in enumerate(muls):
            for cons, deg in muls[ind+1:]:
                if degree == deg and degree not in degrees:
                    constants.append(constant+cons)
                    degrees.append(degree)
                elif degree == deg and degree in degrees:
                    ind = degrees.index(degree)
                    constants[ind] += constant

        for constant, degree in muls:
            if degree not in degrees:
                constants.append(constant)
                degrees.append(degree)

        return Polynomial(constants, degrees)


#TODO FINISH POLYNOMIAL - add division

'''
b = Polynomial([2, 3, 5], [5, 4])
c = Polynomial([5, 5], [5, 4])
x = [x for x in range(-100, 100)]
print(b)
print(c)
print(b*c)
#print(c-b)'''

a = Quadratic(2, 3, 4)
b = a.derivative()
print(b)


