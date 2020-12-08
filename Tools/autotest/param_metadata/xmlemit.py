#!/usr/bin/env python

from xml.sax.saxutils import escape, quoteattr

from emit import Emit
from param import known_param_fields, known_units

INDENT2 = "  "
INDENT4 = "    "
INDENT6 = "      "
INDENT8 = "        "
INDENT10 = "          "

# Emit APM documentation in an machine readable XML format
class XmlEmit(Emit):
    def __init__(self):
        Emit.__init__(self)
        self.wiki_fname = 'apm.pdef.xml'
        self.f = open(self.wiki_fname, mode='w')
        self.preamble = '''<?xml version="1.0" encoding="utf-8"?>
<!-- Dynamically generated list of documented parameters (generated by param_parse.py) -->
'''
        self.f.write(self.preamble)
        self.f.write('''<paramfile>
  <vehicles>\n''')

    def close(self):
        self.f.write(INDENT2 + '</libraries>\n')
        self.f.write('''</paramfile>\n''')
        self.f.close()

    def emit_comment(self, s):
        self.f.write("<!-- " + s + " -->")

    def start_libraries(self):
        self.f.write(INDENT2 + '</vehicles>\n')
        self.f.write(INDENT2 + '<libraries>\n')

    def emit(self, g):
        t = INDENT4 + '''<parameters name=%s>\n''' % quoteattr(g.name)  # i.e. ArduPlane

        for param in g.params:
            # Begin our parameter node
            if hasattr(param, 'DisplayName'):
                t += INDENT6 + '<param humanName=%s name=%s' % (quoteattr(param.DisplayName), quoteattr(param.name))  # i.e. ArduPlane (ArduPlane:FOOPARM)
            else:
                t += INDENT6 + '<param name=%s' % quoteattr(param.name)

            if hasattr(param, 'Description'):
                t += ' documentation=%s' % quoteattr(param.Description)  # i.e. parameter docs
            if hasattr(param, 'User'):
                t += ' user=%s' % quoteattr(param.User)  # i.e. Standard or Advanced

            if hasattr(param, 'Calibration'):
                t += ' calibration=%s' % quoteattr(param.Calibration)  # i.e. Standard or Advanced

            t += ">\n"

            # Add values as chidren of this node
            for field in param.__dict__.keys():
                if field not in ['name', 'DisplayName', 'Description', 'User'] and field in known_param_fields:
                    if field == 'Values' and Emit.prog_values_field.match(param.__dict__[field]):
                        t += INDENT8 + "<values>\n"

                        values = (param.__dict__[field]).split(',')
                        for value in values:
                            v = value.split(':')
                            if len(v) != 2:
                                raise ValueError("Bad value (%s)" % v)
                            t += INDENT10 + '''<value code=%s>%s</value>\n''' % (quoteattr(v[0]), escape(v[1]))  # i.e. numeric value, string label

                        t += INDENT8 + "</values>\n"
                    elif field == 'Units':
                        abreviated_units = param.__dict__[field]
                        if abreviated_units != '':
                            units = known_units[abreviated_units]   # use the known_units dictionary to convert the abreviated unit into a full textual one
                            t += INDENT8 + '''<field name=%s>%s</field>\n''' % (quoteattr(field), escape(abreviated_units))  # i.e. A/s
                            t += INDENT8 + '''<field name=%s>%s</field>\n''' % (quoteattr('UnitText'), escape(units))        # i.e. ampere per second
                    else:
                        t += INDENT8 + '''<field name=%s>%s</field>\n''' % (quoteattr(field), escape(param.__dict__[field]))  # i.e. Range: 0 10

            t += INDENT6 + '''</param>\n'''
        t += INDENT4 + '''</parameters>\n'''

        # print t
        self.f.write(t)