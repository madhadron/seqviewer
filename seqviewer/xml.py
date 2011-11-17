def tag(name):
    def f(classes, *args):
        if not(isinstance(classes, list)):
            classes = [classes]
        return ("""<%s class="%s">""" % (name, ' '.join(classes))) + \
            ''.join(args) + \
            ("""</%s>""" % (name,))
    return f

div = tag('div')
span = tag('span')

def unit_svg(*args):
    return """<svg preserveAspectRatio="none" viewbox="0 -0.05 1 1.05" version="1.1">""" + \
        ''.join(args) + """</svg>"""

def M(x,y):
    return "M%f,%f" % (x,y)

def L(x,y):
    return "L%f,%f" % (x,y)

def path(d, stroke="black", strokeWidth="0.01", fill="none"):
    dstr = ''.join(d)
    return """<path stroke="%s" strokeWidth="%s" fill="%s" d="%s" />""" % \
        (stroke, strokeWidth, fill, dstr)
        
