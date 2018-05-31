from django import template

register = template.Library()


@register.filter
def keyvalue(d, key):    
    return d.get(key, '')