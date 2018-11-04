from django import template

register = template.Library()


@register.filter
def keyvalue(d, key):    
    return d.get(key, '')


@register.filter(name='split')
def split(sample_id):
    tmp = sample_id.split(":")
    tmp[-1] = tmp[-1].replace(" ","_").replace(".","")
    return tmp[-1]
