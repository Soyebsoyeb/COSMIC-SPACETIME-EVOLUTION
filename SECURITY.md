# Security Policy

This is a numerical research code with no network services, no authentication and no user data. The main risks are unsafe deserialisation or dependency vulnerabilities.

* Configuration files are read with `yaml.safe_load` only.
* Please report suspected vulnerabilities privately to [CONTACT EMAIL PLACEHOLDER] rather than in a public issue, and include a minimal reproduction.
* Supported version: the latest release on the default branch.

Scientific errors (wrong equations, incorrect numbers) are not security issues; please open a normal issue for those.
