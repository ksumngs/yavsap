def yavsap_logo() {

tagline = "(Yet Another Viral Subspecies Analysis Pipeline)"
maxlength = tagline.length()

figlet =
"""\
     __   __ ___     ______    _    ____
     \\ \\ / // \\ \\   / / ___|  / \\  |  _ \\
      \\ V // _ \\ \\ / /\\___ \\ / _ \\ | |_) |
       | |/ ___ \\ V /  ___) / ___ \\|  __/
       |_/_/   \\_\\_/  |____/_/   \\_\\_|\
"""

version = "v${workflow.manifest.version}".center(maxlength)

return \
"""\
${figlet}
${tagline}
${version}\
"""
}
